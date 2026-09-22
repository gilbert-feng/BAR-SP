rm(list = ls())

# ------------------------------------------------------------
# Load custom functions
# ------------------------------------------------------------
source_directory <- function(path, pattern = "\\.[Rr]$", recursive = FALSE) {
  files <- list.files(path = path,
                      pattern = pattern,
                      full.names = TRUE,
                      recursive = recursive)
  for (file in files) {
    source(file, encoding = "UTF-8")
    cat("loading: ", file, "\n")
  }
}

source_directory("functions")

library(RSpectra)
library(readr)
library(nleqslv)
library(optimx)
library(Matrix)
library(ggplot2)

# ------------------------------------------------------------
# Seed and data preparation
# ------------------------------------------------------------
set.seed(123)
APP_file <- "data/"
load(paste0(APP_file, "DataList.rda"))

feature <- dataList$X.data
neighbor_list <- dataList$Wmatrices$knn25$nn
n <- dim(neighbor_list)[1]

adj_mat <- matrix(0, nrow = n, ncol = n)
for (i in 1:n) {
  adj_mat[i, neighbor_list[i, ]] <- 1
}
adj_mat <- pmax(adj_mat, t(adj_mat))
A <- adj_mat

Y <- feature[, 1]
Y <- Y * 100
X <- feature[, 2:ncol(feature)]
X <- scale(X, center = TRUE, scale = TRUE)
X <- as.matrix(X)
X <- cbind(1, X)
n <- nrow(X)
p <- ncol(X)

Y <- as.matrix(Y, ncol = 1)
A <- as.matrix(A)
A_temp_num <- as.numeric(A)
dim(A_temp_num) <- dim(A)
A <- A_temp_num

# ------------------------------------------------------------
# Hyper-parameters (K will be selected by IC)
# ------------------------------------------------------------
iter_max <- 1e4
conv_crit <- 1e-3
r <- 0
K_fold <- 2
sig_level <- 0.05
lambda_seq <- 10^seq(log10(1e-3), log10(1), length.out = 100)

select_num <- 3
total_select_num <- select_num * 3
regular <- 0
cv_sqrt <- 0

# ------------------------------------------------------------
# Information criterion for network dimension K
# IC_n(k) = (||P_hat - P_hat_k||_F^2 + kappa_n * k^2) / n^2
# kappa_n = 2 * h_n^2  (Corollary 3.2)
# ------------------------------------------------------------
ic_network_dim <- function(A,
                           K_max = min(50, nrow(A) - 1),
                           h_n = NULL,
                           kappa_n = NULL) {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2
  diag(A) <- 0
  n <- nrow(A)

  eig <- eigen(A, symmetric = TRUE)
  ord <- order(abs(eig$values), decreasing = TRUE)
  vals <- eig$values[ord]

  if (is.null(kappa_n)) {
    if (is.null(h_n)) {
      # Observable Bernoulli calibration (placeholder).
      # Replace by the exact calibration in Online Appendix S6.1.3 if desired.
      p_bar <- mean(A[upper.tri(A)])
      h_n <- 2 * sqrt(n * p_bar * (1 - p_bar)) + sqrt(log(n))
    }
    kappa_n <- 2 * h_n^2
  }

  K_grid <- 0:K_max
  ic <- numeric(length(K_grid))

  for (j in seq_along(K_grid)) {
    k <- K_grid[j]
    fro2 <- if (k >= n) 0 else sum(vals[(k + 1):n]^2)
    ic[j] <- (fro2 + kappa_n * k^2) / n^2
  }

  K_hat <- K_grid[which.min(ic)]

  list(
    IC = data.frame(K = K_grid, IC = ic),
    K_hat = K_hat,
    kappa_n = kappa_n,
    h_n = h_n
  )
}

ic_res <- ic_network_dim(A, K_max = min(50, n - 1), h = 0.6)
K <- ic_res$K_hat

cat("IC-selected network dimension K =", K, "\n")
cat("kappa_n =", ic_res$kappa_n, "\n")

# Save IC values and plot
ic_df <- ic_res$IC
best_df <- ic_df[ic_df$K == K, , drop = FALSE]

p_ic <- ggplot(ic_df, aes(x = K, y = IC)) +
  geom_line(color = "#2166AC", linewidth = 1) +
  geom_point(color = "#2166AC", size = 1.8) +
  geom_point(
    data = best_df,
    aes(x = K, y = IC),
    color = "red",
    size = 3.2,
    inherit.aes = FALSE
  ) +
  geom_vline(
    xintercept = K,
    linetype = "dashed",
    color = "red",
    linewidth = 0.6
  ) +
  annotate(
    "text",
    x = K,
    y = best_df$IC,
    label = paste0("K = ", K),
    vjust = -1.2,
    hjust = ifelse(K > max(ic_df$K) / 2, 1.1, -0.1),
    color = "red",
    size = 4
  ) +
  scale_x_continuous(breaks = pretty(ic_df$K, n = 8)) +
  labs(
    title = "Information criterion for network dimension",
    x = "Network dimension K",
    y = "IC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

# ------------------------------------------------------------
# Preliminary SVD and heteroskedasticity test
# ------------------------------------------------------------
temp.ProdSVD.fit <- ProdSVD(X, A, K)
ProdSVD.fit <- ProdSVD(X, A, K, r = temp.ProdSVD.fit$r)
temp_X_proj <- ProdSVD.fit$x_proj
temp_Z <- ProdSVD.fit$Z

SP.est <- SP.Prod(X, ProdSVD.fit$x_proj, ProdSVD.fit$Z_hat, Y,
                  ProdSVD.fit$W_cup, K, ProdSVD.fit$r, 1, sig_level)
est_error2 <- SP.est$est_error^2
new_X <- cbind(X, X[, -1]^2)
hettest <- lm(est_error2 ~ new_X - 1)
homo <- ifelse(any(summary(hettest)$coefficients[, 4] < 0.1), 0, 1)

# ------------------------------------------------------------
# Create independent result folder (no longer under "results")
# Top-level: IC_K_selection, then homo / hete
# ------------------------------------------------------------
base_dir <- "IC_K_selection"
prefix <- ifelse(homo == 1, paste0(base_dir, "/homo/"), paste0(base_dir, "/hete/"))
dir_name <- paste0(prefix, "n = ", n, ", K = ", K, "/")
dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)

# Save IC plot and IC values into the final folder
ggsave(
  filename = paste0(dir_name, "IC_K_selection.pdf"),
  plot = p_ic,
  device = "pdf",
  width = 7,
  height = 5
)
write.csv(ic_res$IC, paste0(dir_name, "IC_K_selection.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Initialize estimation results matrix
# ------------------------------------------------------------
est_metric <- matrix(0, 6, 4)

# (1) BAR-SP method
BAR.SP.est <- BAR.SP.Prod(
  X, ProdSVD.fit$x_proj, ProdSVD.fit$Z_hat, Y, ProdSVD.fit$W_cup,
  K, ProdSVD.fit$r, homo, K_fold, lambda_seq, sig_level,
  regular = regular, cv_sqrt = cv_sqrt
)
active_ratio <- sum(abs(BAR.SP.est$beta) > 1e-5) / p
est_metric[, 1] <- c(
  BAR.SP.est$chisq.p, BAR.SP.est$r, active_ratio,
  BAR.SP.est$lambda_opt, BAR.SP.est$out_mse,
  BAR.SP.est$in_mse
)
alpha_hat <- BAR.SP.est$alpha

temp_beta <- BAR.SP.est$beta
top_inactive <- head(order(abs(temp_beta), decreasing = FALSE), select_num)
top_active <- which(abs(BAR.SP.est$beta) > 1e-5)
top_indices <- c(top_active, top_inactive)
param_sig_num <- matrix(-10, length(top_indices) + 1, 12)

param_sig_num[, 1] <- c(-10, BAR.SP.est$beta[top_indices])
active_colnames <- paste0("V", top_active)
temp_p <- rep(-10, length(top_indices))
for (i in 1:length(top_active)) {
  temp_p[i] <- BAR.SP.est$act.coef.mat[
    rownames(BAR.SP.est$act.coef.mat) == active_colnames[i], 4
  ]
}
param_sig_num[, 2] <- c(BAR.SP.est$chisq.p, temp_p)
param_sig_num[, 3] <- c(-10, top_indices)

# (2) SP method
SP.est <- SP.Prod.pred(
  X, ProdSVD.fit$x_proj, ProdSVD.fit$Z_hat, Y, ProdSVD.fit$W_cup,
  K, ProdSVD.fit$r, 1, K_fold, sig_level, regular = regular
)
active_ratio <- sum(abs(SP.est$beta) > 1e-5) / p
est_metric[, 2] <- c(
  SP.est$chisq.p, SP.est$r, active_ratio,
  0, SP.est$out_mse,
  SP.est$in_mse
)

param_sig_num[, 4] <- c(-10, SP.est$beta[top_indices])
active_colnames <- paste0("V", top_indices)
for (i in 1:length(top_indices)) {
  temp_p[i] <- SP.est$act.coef.mat[
    rownames(SP.est$act.coef.mat) == active_colnames[i], 4
  ]
}
param_sig_num[, 5] <- c(SP.est$chisq.p, temp_p)
param_sig_num[, 6] <- c(-10, top_indices)

# Compare significant variables between BAR-SP and SP
SP_sig_var <- rownames(SP.est$act.coef.mat[SP.est$act.coef.mat[, 4] < 0.05, ])
BARSP_sig_var <- rownames(BAR.SP.est$act.coef.mat[BAR.SP.est$act.coef.mat[, 4] < 0.05, ])

act_diff_var <- setdiff(SP_sig_var, BARSP_sig_var)
act_diff_var <- parse_number(act_diff_var)
write.csv(colnames(X)[act_diff_var],
          paste0(dir_name, "BAR not sig SP sig.csv"), row.names = FALSE)

act_diff_var <- setdiff(BARSP_sig_var, SP_sig_var)
act_diff_var <- parse_number(act_diff_var)
write.csv(colnames(X)[act_diff_var],
          paste0(dir_name, "BAR sig SP not sig.csv"), row.names = FALSE)

# Group structure investigation
select_beta <- BAR.SP.est$beta[top_active]
select_beta_order <- order(select_beta, decreasing = TRUE)
select_beta_diff <- diff(select_beta[select_beta_order])
group_dis <- 0.01
select_ind <- which((abs(select_beta_diff) < group_dis) & (select_beta_diff != 0))

temp_name <- matrix(0, 2, length(select_ind))
temp_name[1, ] <- colnames(X)[top_indices[select_beta_order[select_ind]]]
temp_name[2, ] <- colnames(X)[top_indices[select_beta_order[select_ind + 1]]]
temp_name <- t(temp_name)
rownames(temp_name) <- select_ind
write.csv(temp_name, paste0(dir_name, "group effect.csv"), row.names = TRUE)

# Network effect terms
show_num <- 5
beta_thres <- 0
pos_alpha_ind <- order(alpha_hat, decreasing = TRUE)[1:show_num]
neg_alpha_ind <- order(alpha_hat, decreasing = FALSE)[1:show_num]
pos_coef_ind <- which(BAR.SP.est$beta > beta_thres)
neg_coef_ind <- which(BAR.SP.est$beta < -beta_thres)

agg_var_names <- colSums(X[pos_alpha_ind, ])
agg_ind <- which(agg_var_names > 0)
write.csv(colnames(X)[agg_ind],
          paste0(dir_name, "high network effect terms.csv"), row.names = TRUE)

agg_inter_ind <- intersect(agg_ind, pos_coef_ind)
write.csv(colnames(X)[agg_inter_ind],
          paste0(dir_name, "high network effect and high coef terms.csv"), row.names = TRUE)

agg_inter_ind <- intersect(agg_ind, neg_coef_ind)
write.csv(colnames(X)[agg_inter_ind],
          paste0(dir_name, "high network effect and low coef terms.csv"), row.names = TRUE)

agg_var_names <- colSums(X[neg_alpha_ind, ])
agg_ind <- which(agg_var_names > 0)
write.csv(colnames(X)[agg_ind],
          paste0(dir_name, "low network effect terms.csv"), row.names = TRUE)

agg_inter_ind <- intersect(agg_ind, pos_coef_ind)
write.csv(colnames(X)[agg_inter_ind],
          paste0(dir_name, "low network effect and high coef terms.csv"), row.names = TRUE)

agg_inter_ind <- intersect(agg_ind, neg_coef_ind)
write.csv(colnames(X)[agg_inter_ind],
          paste0(dir_name, "low network effect and low coef terms.csv"), row.names = TRUE)

# (3) OLS
OLS.est <- OLSCV(X, Y, homo)
active_ratio <- sum(abs(OLS.est$beta_hat) > 1e-5) / p
est_metric[, 3] <- c(
  0, 0, active_ratio,
  0, 0,
  OLS.est$in_mse
)

param_sig_num[, 7] <- c(-10, OLS.est$beta_hat[top_indices])
param_sig_num[, 8] <- c(-10, 1 - pnorm(abs(OLS.est$beta_hat / diag(OLS.est$beta_sd))[top_indices]))
param_sig_num[, 9] <- c(-10, top_indices)

# (4) SIM
SIM.est <- SIMCV(X, Y, A, homo)
active_ratio <- sum(abs(SIM.est$phi_est[-1]) > 1e-5) / p
SIM.p.values <- 1 - pnorm(abs(SIM.est$phi_est / SIM.est$phi_sd))
est_metric[, 4] <- c(
  SIM.p.values[1], 0, active_ratio,
  0, 0,
  SIM.est$in_mse
)
param_sig_num[, 10] <- SIM.est$phi_est[c(0, top_indices) + 1]
param_sig_num[, 11] <- SIM.p.values[c(0, top_indices) + 1]
param_sig_num[, 12] <- c(-10, top_indices)

# Final outputs
write.csv(est_metric, paste0(dir_name, "metric.csv"))
feature_name <- colnames(X)
results_time <- data.frame(
  param_sig_num,
  feature_name = c("NetEff", feature_name[top_indices])
)
write.csv(results_time, paste0(dir_name, "param10.csv"), row.names = FALSE)
write.csv(BAR.SP.est$project_measures,
          paste0(dir_name, "proj measure.csv"), row.names = FALSE)

# Alpha plot data
sorted_alpha <- sort(alpha_hat, decreasing = TRUE)
n_alpha_ind <- order(alpha_hat, decreasing = TRUE)[(n - 9):n]

alpha_df <- data.frame(
  index = 1:n,
  alpha = alpha_hat,
  group = ifelse(1:n %in% n_alpha_ind, "Top 10", "Other")
)
write.csv(alpha_df, paste0(dir_name, "alpha.csv"), row.names = FALSE)

xi_df <- data.frame(
  index = 1:n,
  alpha = BAR.SP.est$xi_hat,
  group = ifelse(1:n %in% n_alpha_ind, "Top 10", "Other")
)
write.csv(BAR.SP.est$xi_hat, paste0(dir_name, "xi.csv"), row.names = FALSE)
write.csv(BAR.SP.est$fit.residual,
          paste0(dir_name, "BAR residual.csv"), row.names = FALSE)

cat("All results saved in:", dir_name, "\n")
