rm(list=ls())
source_directory <- function(path, pattern = "\\.[Rr]$", recursive = FALSE) {
  files <- list.files(path = path,
                      pattern = pattern,
                      full.names = TRUE,
                      recursive = recursive)
  for (file in files) {
    source(file, encoding = 'UTF-8')
    cat("loading: ", file, "\n") 
  }
}

source_directory("functions")
library(RSpectra)
library(readr)
library(nleqslv)
library(optimx)
library(Matrix)

# Set seed
set.seed(123)
APP_file <- "data/"
load(paste0(APP_file,"DataList.rda"))
feature <- dataList$X.data
neighbor_list <- dataList$Wmatrices$knn25$nn
n <- dim(neighbor_list)[1]
adj_mat <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  adj_mat[i, neighbor_list[i,]] <- 1
}

# Make symmetric
adj_mat <- pmax(adj_mat, t(adj_mat))
A <- adj_mat

Y <- feature[,1]
Y <- Y*100
X <- feature[,2:ncol(feature)]
X <- scale(X, center = TRUE, scale = TRUE)

X <- as.matrix(X)

# Augment covariates with intercept
X <- cbind(1,X)
n <- nrow(X)
p <- ncol(X)


Y <- as.matrix(Y,ncol = 1)
A <- as.matrix(A)
A_temp_num <- as.numeric(A)
dim(A_temp_num) <- dim(A)
A <- A_temp_num


# Define hyperparameters
K <- 25
iter_max <- 1e4
conv_crit <- 1e-3
r <- 0  
K_fold <- 2
sig_level <- 0.05
lambda_seq <- NULL
lambda_seq <- 10^seq(log10(1e-3), log10(1), length.out = 100)

select_num <- 3
total_select_num <- select_num*3
regular <- 0
cv_sqrt <- 0


#############################
# 1. alpha cluster analysis
#############################
library(ggplot2)
library(dplyr)
library(maps)
library(ggrepel)   # For smart label avoidance to ensure full name display

# ==================== 1. Prepare data ====================
set.seed(123)
results_prefix <- "results/hete"
alpha_hat <- read.csv(paste0(results_prefix, "/alpha.csv"))
paper_id <- 1:n

# Assume Y is average growth (ensure Y is defined)
data <- data.frame(
  id = 1:n,
  x = Y,                       # Average growth
  y = alpha_hat$alpha          # Network effect
)

# City coordinates (n×2)
city_location <- dataList$Wmatrices$knn25$x
colnames(city_location) <- c("lon", "lat")

# Whether it is a capital (feature$Capital should be 0/1)
data$is_capital <- feature$Capital
data$cluster <- ifelse(data$is_capital == 1, "Capital", "Non-Capital")

# Adjacency matrix A
edges <- as.data.frame(which(A != 0, arr.ind = TRUE, useNames = FALSE))
colnames(edges) <- c("from", "to")

# ==================== 2. Identify top five capitals (descending by average growth) ====================
capital_ids <- data$id[data$is_capital == 1]
capital_data <- data[capital_ids, ]
top5_capitals <- capital_data %>%
  arrange(desc(x)) %>%        # Sort descending by average growth
  slice_head(n = 5) %>%
  pull(id)

# ==================== 3. Filter edges (only those involving top five capitals) ====================
edges_filtered <- edges %>%
  filter(from %in% top5_capitals | to %in% top5_capitals)

# Prepare coordinates for the edges (scatter plot coordinates)
plot_links <- edges_filtered %>%
  left_join(data, by = c("from" = "id")) %>%
  rename(x1 = x, y1 = y) %>%
  left_join(data, by = c("to" = "id")) %>%
  rename(x2 = x, y2 = y)

# ==================== 4. Match country names for top five capitals (based on coordinates) ====================
# Function: return country name from coordinates (consistent with maps package region)
get_country_name <- function(lon, lat) {
  country <- map.where("world", lon, lat)
  if (is.na(country)) return(NA)
  # Remove ":Alaska" part from "USA:Alaska"
  country <- strsplit(country, ":")[[1]][1]
  return(country)
}

# Extract coordinates of top five capitals
top5_coords <- city_location[top5_capitals, , drop = FALSE]
top5_countries <- apply(top5_coords, 1, function(row) {
  get_country_name(row[1], row[2])
})

# Add country names to data frame (only for top five capitals)
data$country <- NA
data$country[top5_capitals] <- top5_countries

# ==================== 5. Draw scatter plot ====================
p <- ggplot() +
  # Edges (only those connected to top five capitals)
  geom_segment(data = plot_links,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "gray50", alpha = 0.3, linewidth = 0.3) +
  # Scatter points: color/shape by Capital / Non-Capital
  geom_point(data = data,
             aes(x = x, y = y, color = cluster, shape = cluster),
             size = 2, alpha = 0.9) +
  # Add country name labels for top five capitals (using ggrepel for auto avoidance + margin adjustment)
  geom_text_repel(data = data %>% filter(id %in% top5_capitals),
                  aes(x = x, y = y, label = country),
                  size = 3.5, color = "black", fontface = "bold",
                  box.padding = 0.5, point.padding = 0.3,
                  force = 2, seed = 123, show.legend = FALSE) +
  labs(x = "Average growth", y = "Network effect",
       color = "Type", shape = "Type") +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        # Increase plot margins to prevent labels from being clipped
        plot.margin = margin(15, 15, 15, 15))

# Save plot: increase width/height and disable clipping to ensure full text display
ggsave("stat_alpha_cluster_top5.pdf", 
       plot = p, 
       device = cairo_pdf, 
       width = 10,      # increased from 8 to 10
       height = 7,      # increased from 6 to 7
       units = "in")



#############################
# 2. Disturbance analysis
#############################

disturb_df <- read.csv("disturb agg.csv")
#disturb_df <- disturb_df[seq(2,10,2),]
row_df <- nrow(disturb_df)
col_names <- c("BAR-SP", "SP", "SIM", "0.05 level")
df <- tibble(
  x = rep(disturb_df$degree, length(col_names)),
  y = c(disturb_df[,2], disturb_df[,3], disturb_df[,4], rep(0.05,row_df)),
  Group = factor(rep(col_names, each = row_df),
                 levels = col_names)
)

p <- ggplot(df, aes(x = x, y = y, color = Group, linetype = Group, shape = Group)) +
  # Add lines
  geom_line(linewidth = 0.8) +
  
  # Add points
  geom_point(size = 3, fill = "white") +
  
  # Axis labels and title
  labs(
    x = "Degree of disturbance",
    y = ""
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",         # Place legend above the plot
    legend.direction = "horizontal", # Arrange legend items horizontally
    legend.key.width = unit(1.5, "cm"),  
    # Optional: fine-tune spacing between legend text and graphics
    legend.spacing.x = unit(0.2, "cm")
  )

ggsave("stat disturb.pdf", 
       plot = p,
       device = cairo_pdf,  # Use cairo device
       width = 8, 
       height = 6,
       units = "in")

#############################
# 3. Alpha geological analysis
#############################
# Load necessary packages
library(ggplot2)
library(maps)
library(viridis)
library(dplyr)
library(ggforce)   # <--- Change: load ggforce for drawing circles


# Construct data (according to your actual structure)
df <- cbind(alpha_hat$alpha, dataList$Wmatrices$knn25$x)
colnames(df) <- c("network", "lon", "lat")

# Check data range
summary(df[, c("lon", "lat")])

# ------------------------------
# Get European basemap
# ------------------------------
world <- map_data("world")
europe <- world %>%
  filter(long >= -25 & long <= 45 & lat >= 35 & lat <= 70) %>%
  filter(!(region %in% c("Russia")))

# Compute centroid for each country (for labeling)
country_labels <- europe %>%
  group_by(region) %>%
  summarise(
    long = mean(range(long)),   # centroid longitude (could also use median or mean)
    lat  = mean(range(lat)),    # centroid latitude
    .groups = "drop"
  ) %>%
  # Optional: filter out countries that are too small or not in main display area
  filter(!region %in% c("Isle of Man", "Faeroe Islands", "Svalbard"))

# ------------------------------
# Draw spatial distribution plot
# ------------------------------
p <- ggplot() +
  # European basemap
  geom_polygon(data = europe,
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "black", linewidth = 0.2) +
  # Add country name labels (using centroid coordinates)
  ggrepel::geom_text_repel(data = country_labels,
                           aes(x = long, y = lat, label = region),
                           size = 2.5, colour = "black", fontface = "bold") +
  # Overlay data points, color mapped to network effect
  geom_point(data = df,
             aes(x = lon, y = lat, color = network),
             size = 2, alpha = 0.8) +
  # Use viridis palette and reverse direction: larger alpha -> darker color
  scale_color_viridis_c(option = "plasma", direction = -1, name = expression(network)) +
  # <--- Change start: add two circle markers ------------------
# Circle 1: mark area around Ireland (center ≈ (-8, 53), radius 1.5 degrees)
geom_circle(aes(x0 = -8, y0 = 53, r = 1.5),
            color = "red", size = 0.8, linetype = "dashed", fill = NA) +
  # Circle 2: mark area around Estonia, Latvia, Lithuania, Poland (center ≈ (20, 56), radius 4 degrees)
  geom_circle(aes(x0 = 22, y0 = 56, r = 5),
              color = "blue", size = 0.8, linetype = "dashed", fill = NA) +
  # <--- Change end ----------------------------------
# Set coordinate limits (focus on Europe)
coord_quickmap(xlim = c(-12, 32), ylim = c(35, 66)) +
  # Titles and theme
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

# Save plot
ggsave("alpha dist.pdf", 
       plot = p,
       device = cairo_pdf,
       width = 8, 
       height = 6,
       units = "in")

#############################
# 4. Residual geological analysis
#############################
residual <- read.csv(paste0(results_prefix,"/BAR residual.csv"))

df <- cbind(residual, dataList$Wmatrices$knn25$x)
colnames(df) <- c("residual", "lon", "lat")

p <- ggplot() +
  # European basemap
  geom_polygon(data = europe,
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "black", linewidth = 0.2) +
  # Add country name labels (using centroid coordinates)
  ggrepel::geom_text_repel(data = country_labels,
                           aes(x = long, y = lat, label = region),
                           size = 2.5, colour = "black", fontface = "bold") +
  # Overlay data points, color mapped to residual
  geom_point(data = df,
             aes(x = lon, y = lat, color = residual),
             size = 2, alpha = 0.8) +
  # Use viridis palette and reverse direction: larger residual -> darker color
  scale_color_viridis_c(option = "viridis", direction = -1, name = expression(residual)) +
  # <--- Change start: add circle markers ------------------
# Circle 1: mark area around Ireland (center ≈ (-8, 53), radius 1.5 degrees)
geom_circle(aes(x0 = -8, y0 = 53, r = 1.5),
            color = "red", size = 0.8, linetype = "dashed", fill = NA) +
  # <--- Change end ----------------------------------
# Set coordinate limits (focus on Europe)
coord_quickmap(xlim = c(-12, 32), ylim = c(35, 66)) +
  # Titles and theme
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

# Save plot
ggsave("residual dist.pdf", 
       plot = p,
       device = cairo_pdf,
       width = 8, 
       height = 6,
       units = "in")

#############################
# 5. Xi geological analysis
#############################
residual <- read.csv(paste0(results_prefix,"/xi.csv"))

df <- cbind(residual, dataList$Wmatrices$knn25$x)
colnames(df) <- c("common", "lon", "lat")

p <- ggplot() +
  # European basemap
  geom_polygon(data = europe,
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "black", linewidth = 0.2) +
  # Add country name labels (using centroid coordinates)
  ggrepel::geom_text_repel(data = country_labels,
                           aes(x = long, y = lat, label = region),
                           size = 2.5, colour = "black", fontface = "bold") +
  # Overlay data points, color mapped to common component
  geom_point(data = df,
             aes(x = lon, y = lat, color = common),
             size = 2, alpha = 0.8) +
  # Use viridis palette and reverse direction: larger common -> darker color
  scale_color_viridis_c(option = "viridis", direction = -1, name = expression(common)) +
  # <--- Change start: add two circle markers ------------------
# Circle 1: mark area around Ireland (center ≈ (-8, 53), radius 1.5 degrees)
geom_circle(aes(x0 = -8, y0 = 53, r = 1.5),
            color = "red", size = 0.8, linetype = "dashed", fill = NA) +
  # Circle 2: mark area around Estonia, Latvia, Lithuania, Poland (center ≈ (20, 56), radius 4 degrees)
  geom_circle(aes(x0 = 22, y0 = 56, r = 5),
              color = "blue", size = 0.8, linetype = "dashed", fill = NA) +
  # <--- Change end ----------------------------------
# Set coordinate limits (focus on Europe)
coord_quickmap(xlim = c(-12, 32), ylim = c(35, 66)) +
  # Titles and theme
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")

# Save plot
ggsave("xi dist.pdf", 
       plot = p,
       device = cairo_pdf,
       width = 8, 
       height = 6,
       units = "in")

