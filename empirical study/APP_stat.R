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

# Seed
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

# symmetric
adj_mat <- pmax(adj_mat, t(adj_mat))
A <- adj_mat

Y <- feature[,1]
Y <- Y*100
X <- feature[,2:ncol(feature)]
X <- scale(X, center = TRUE, scale = TRUE)

X <- as.matrix(X)

# covariates augmented
X <- cbind(1,X)
n <- nrow(X)
p <- ncol(X)


Y <- as.matrix(Y,ncol = 1)
A <- as.matrix(A)
A_temp_num <- as.numeric(A)
dim(A_temp_num) <- dim(A)
A <- A_temp_num


# define hyper-parameters
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



temp.ProdSVD.fit <- ProdSVD(X,A,K)
# ProdSVD.my.fit <- my_boot(A, X, K, 50)
# ProdSVD.my.fit$r
# ProdSVD.my.fit$thr
temp.ProdSVD.fit$r

# preliminary SVD 
ProdSVD.fit <- ProdSVD(X,A,K,r=temp.ProdSVD.fit$r) # r=6 is estimated with bootstrap method
temp_X_proj <- ProdSVD.fit$x_proj
temp_Z <- ProdSVD.fit$Z

# a simple method to determine the existence of heteroskedastic disturbance
SP.est <- SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,1,sig_level)
est_error2 <- SP.est$est_error^2
new_X <- cbind(X, X[, -1]^2)  # add squared terms
hettest <- lm(est_error2 ~ new_X-1)
homo <- ifelse(any(summary(hettest)$coefficients[, 4] < 0.1), 0, 1)

#homo <- 1

# initialize estimation results matrix
est_metric <- matrix(0, 6, 4)

prefix <- ifelse(homo == 1, "results/homo/", "results/hete/")
dir_name <- paste0(prefix, "n = ", n)
dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)


# (1) A-SP method
BAR.SP.est <- BAR.SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,homo,K_fold,lambda_seq,sig_level,regular=regular,cv_sqrt)
active_ratio <- sum(abs(BAR.SP.est$beta)>1e-5)/p
est_metric[, 1] <- c(BAR.SP.est$chisq.p, BAR.SP.est$r, active_ratio, 
                     BAR.SP.est$lambda_opt, BAR.SP.est$out_mse, 
                     BAR.SP.est$in_mse)
alpha_hat <- BAR.SP.est$alpha

# select significant features (without intercept)
temp_beta <- BAR.SP.est$beta
#temp_beta[1] <- NA
# active_pos_beta <- temp_beta
# active_pos_beta[active_pos_beta<0] <- 0
# top_active_pos <- head(order(active_pos_beta, decreasing = TRUE), select_num)
# active_neg_beta <- temp_beta
# active_neg_beta[active_neg_beta>0] <- 0
# top_active_neg <- head(order(active_neg_beta, decreasing = FALSE), select_num)
top_inactive <- head(order(abs(temp_beta), decreasing = FALSE), select_num)
top_active <- which(abs(BAR.SP.est$beta)>1e-5)
top_indices <- c(top_active, top_inactive)
param_sig_num <- matrix(-10, length(top_indices)+1, 12)


param_sig_num[, 1] <- c(-10,BAR.SP.est$beta[top_indices])
active_colnames <- paste0("V",top_active)
temp_p <- rep(-10,length(top_indices))
for (i in 1:length(top_active)) {
  temp_p[i] <- BAR.SP.est$act.coef.mat[rownames(BAR.SP.est$act.coef.mat)==active_colnames[i],4]
}
param_sig_num[, 2] <- c(BAR.SP.est$chisq.p, temp_p)
param_sig_num[, 3] <- c(-10, top_indices)

# (2) SP 
SP.est <- SP.Prod.pred(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,1,K_fold,sig_level,regular=regular)
active_ratio <- sum(abs(SP.est$beta)>1e-5)/p
est_metric[, 2] <- c(SP.est$chisq.p, SP.est$r, active_ratio, 
                     0, SP.est$out_mse, 
                     SP.est$in_mse)

param_sig_num[, 4] <- c(-10, SP.est$beta[top_indices])
active_colnames <- paste0("V",top_indices)
for (i in 1:length(top_indices)) {
  temp_p[i] <- SP.est$act.coef.mat[rownames(SP.est$act.coef.mat)==active_colnames[i],4]
}
param_sig_num[, 5] <- c(SP.est$chisq.p, temp_p)
param_sig_num[, 6] <- c(-10, top_indices)

## compare the results between BAR-SP and SP
### significant varibles comparison
SP_sig_var <- rownames(SP.est$act.coef.mat[SP.est$act.coef.mat[,4]<0.05,])
BARSP_sig_var <- rownames(BAR.SP.est$act.coef.mat[BAR.SP.est$act.coef.mat[,4]<0.05,])
length(SP_sig_var)
length(BARSP_sig_var)

act_diff_var <- setdiff(SP_sig_var,BARSP_sig_var)
act_diff_var <- parse_number(act_diff_var)
write.csv(colnames(X)[act_diff_var], paste0(prefix, "BAR not sig SP sig.csv"), row.names = F)


act_diff_var <- setdiff(BARSP_sig_var,SP_sig_var)
act_diff_var <- parse_number(act_diff_var)
write.csv(colnames(X)[act_diff_var], paste0(prefix, "BAR sig SP not sig.csv"), row.names = F)

### coefficient group structure investigation
select_beta <- BAR.SP.est$beta[top_active]
select_beta_order <- order(select_beta,decreasing = T)
select_beta_diff <- diff(select_beta[select_beta_order])
group_dis <- 0.01
select_ind <- which((abs(select_beta_diff)<group_dis)&(select_beta_diff!=0))

temp_name <- matrix(0,2,length(select_ind))
temp_name[1,] <- colnames(X)[top_indices[select_beta_order[select_ind]]]
temp_name[2,] <- colnames(X)[top_indices[select_beta_order[select_ind+1]]]
temp_name <- t(temp_name)
rownames(temp_name) <- select_ind
write.csv(temp_name, paste0(prefix, "group effect.csv"), row.names = T)


# coef sig + alpha
show_num <- 5
beta_thres <- 0
pos_alpha_ind <- order(alpha_hat, decreasing = TRUE)[1:show_num]
neg_alpha_ind <- order(alpha_hat, decreasing = FALSE)[1:show_num]
# pos_coef_ind <- order(BAR.SP.est$beta, decreasing = TRUE)[1:show_num]
# neg_coef_ind <- order(BAR.SP.est$beta, decreasing = FALSE)[1:show_num]
pos_coef_ind <- which(BAR.SP.est$beta>beta_thres)
neg_coef_ind <- which(BAR.SP.est$beta< -beta_thres)


## high network effect terms
agg_var_names <- colSums(X[pos_alpha_ind,])
agg_ind <- which(agg_var_names>0)
write.csv(colnames(X)[agg_ind], paste0(prefix, "high network effect terms.csv"), row.names = T)

## high network effect and high coef terms 
agg_inter_ind <- intersect(agg_ind,pos_coef_ind)
write.csv(colnames(X)[agg_inter_ind], paste0(prefix, "high network effect and high coef terms.csv"), row.names = T)

## high network effect and low coef terms 
agg_inter_ind <- intersect(agg_ind,neg_coef_ind)
write.csv(colnames(X)[agg_inter_ind], paste0(prefix, "high network effect and low coef terms.csv"), row.names = T)


## low network effect terms
agg_var_names <- colSums(X[neg_alpha_ind,])
agg_ind <- which(agg_var_names>0)
write.csv(colnames(X)[agg_ind], paste0(prefix, "low network effect terms.csv"), row.names = T)

## low network effect and high coef terms 
agg_inter_ind <- intersect(agg_ind,pos_coef_ind)
write.csv(colnames(X)[agg_inter_ind], paste0(prefix, "low network effect and high coef terms.csv"), row.names = T)

## high network effect and low coef terms 
agg_inter_ind <- intersect(agg_ind,neg_coef_ind)
write.csv(colnames(X)[agg_inter_ind], paste0(prefix, "low network effect and low coef terms.csv"), row.names = T)



# logistic -> causal
# criterion -> forecast
# 
# colnames(X)[param_sig_num[,3]]
# 
# 
# paper_id$title[which(alpha_hat>6)]
# 
# 
# paper_id$title[which(alpha_hat<(-4.5))] # diversify; no pattern

# (3) OLS
OLS.est <- OLSCV(X,Y,homo)
active_ratio <- sum(abs(OLS.est$beta_hat)>1e-5)/p
est_metric[, 3] <- c(0, 0, active_ratio, 
                     0, 0, 
                     OLS.est$in_mse)

param_sig_num[, 7] <- c(-10, OLS.est$beta_hat[top_indices])
param_sig_num[, 8] <- c(-10, 1-pnorm(abs(OLS.est$beta_hat / diag(OLS.est$beta_sd))[top_indices]))
param_sig_num[, 9] <- c(-10, top_indices)

# (4) SIM
# row_sums <- rowSums(A)
# A_row_norm <- sweep(A, 1, row_sums, FUN = "/")
SIM.est <- SIMCV(X,Y,A,homo)
active_ratio <- sum(abs(SIM.est$phi_est[-1])>1e-5)/p
SIM.p.values <- 1-pnorm(abs(SIM.est$phi_est / SIM.est$phi_sd))
est_metric[, 4] <- c(SIM.p.values[1], 0, active_ratio, 
                     0, 0, 
                     SIM.est$in_mse)
param_sig_num[, 10] <- SIM.est$phi_est[c(0,top_indices)+1] # autoregressive term
param_sig_num[, 11] <- SIM.p.values[c(0,top_indices)+1] # autoregressive term
param_sig_num[, 12] <- c(-10,top_indices)


# final output
write.csv(est_metric, paste0(prefix, "metric.csv"))
feature_name <- colnames(X)
results_time <- data.frame(
  param_sig_num,
  feature_name = c("NetEff", feature_name[top_indices])
)
write.csv(results_time, paste0(prefix, "param10.csv"), row.names = FALSE)
write.csv(BAR.SP.est$project_measures, paste0(prefix, "proj measure.csv"), row.names = FALSE)


# alpha plot 
library(ggplot2)
sorted_alpha <- sort(alpha_hat, decreasing = TRUE)
n_alpha_ind <- order(alpha_hat, decreasing = TRUE)[(n-9):n]

alpha_df <- data.frame(
  index = 1:n,
  alpha = alpha_hat,
  group = ifelse(1:n %in% n_alpha_ind, "Top 10", "Other")
)

# ggplot(alpha_df, aes(x = index, y = alpha, color = group, shape = group)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = c("Top 10" = "red", "Other" = "blue")) +
#   scale_shape_manual(values = c("Top 10" = 19, "Other" = 8)) +
#   labs(y = "Alpha", title = "Estimated Network Effect") +
#   theme_minimal() +
#   theme(legend.title = element_blank())

write.csv(alpha_df, paste0(prefix, "alpha.csv"), row.names = FALSE)

xi_df <- data.frame(
  index = 1:n,
  alpha = BAR.SP.est$xi_hat,
  group = ifelse(1:n %in% n_alpha_ind, "Top 10", "Other")
)
write.csv(BAR.SP.est$xi_hat, paste0(prefix, "xi.csv"), row.names = FALSE)


write.csv(BAR.SP.est$fit.residual, paste0(prefix, "BAR residual.csv"), row.names = FALSE)


# network perturbation experiment
min_degree <- 0
max_degree <- 0.1 # 0.0005
disturb_degree <- seq(min_degree, max_degree, length.out = 10)
n_rep <- 50
results <- matrix(0, length(disturb_degree), 3)
results[1, ] <- c(
  est_metric[c(1,3), 1],
  sum(BAR.SP.est$act.coef.mat[,4]<sig_level)/p)

results_SP <- matrix(0, length(disturb_degree), 3)
results_SP[1, ] <- c(
  est_metric[c(1,3), 2],
  sum(SP.est$act.coef.mat[,4]<sig_level)/p)

results_SIM <- matrix(0, length(disturb_degree), 3)
results_SIM[1, ] <- c(
  est_metric[c(1,3), 4],
  sum(SIM.p.values[-1]<sig_level)/p)


disturb_network <- function(P, degree) {
  A <- P
  mask <- upper.tri(A)
  upper_part <- A[mask]
  
  idx_1 <- (upper_part == 1)
  rand_mask_0 <- runif(length(upper_part)) <= degree
  upper_part[idx_1 & rand_mask_0] <- 0
  
  idx_0 <- (upper_part == 0)
  rand_mask_1 <- runif(length(upper_part)) <= degree
  upper_part[idx_0 & rand_mask_1] <- 1
  
  A[mask] <- upper_part
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  return(A)
}


for (i in seq_along(disturb_degree)) {
  degree <- disturb_degree[i]

  if (i >= 2) {
    temp_results <- matrix(0, n_rep, 3)
    temp_results_SP <- matrix(0, n_rep, 3)
    temp_results_SIM <- matrix(0, n_rep, 3)

    for (j in 1:n_rep) {
      A_temp <- disturb_network(A, degree)
      eigen_A <- eigs_sym(A=A_temp,k=K)
      W <- matrix(eigen_A$vectors[,1:K],ncol=K)
      temp_ProdSVD <- svd(t(temp_Z)%*%W,nu=p,nv=K)
      S_hat <- temp_ProdSVD$d
      Z_hat <- temp_Z%*%temp_ProdSVD$u
      W_cup <- W%*%temp_ProdSVD$v
      ProdSVD.fit <- list(x_proj = temp_X_proj, S_hat = S_hat,
                          Z = temp_Z, Z_hat = Z_hat,
                          W = W, W_cup = W_cup,
                          r = r)
      BAR.SP.est <- BAR.SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,homo,K_fold,lambda_seq,sig_level,regular=regular,cv_sqrt = cv_sqrt)
      SP.est <- SP.Prod.pred(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,1,K_fold,sig_level,regular=regular)
      SIM.est <- SIMCV(X,Y,A_temp,homo)
      
      active_ratio <- sum(abs(SIM.est$phi_est[-1])>1e-5)/p
      SIM.p.values <- 1-pnorm(abs(SIM.est$phi_est / SIM.est$phi_sd))

      
      temp_results[j, ] <- c(
        BAR.SP.est$chisq.p,
        sum(abs(BAR.SP.est$beta)>1e-5)/p,
        sum(BAR.SP.est$act.coef.mat[,4]<sig_level)/p)
      temp_results_SP[j, ] <- c(
        SP.est$chisq.p,
        sum(abs(SP.est$beta)>1e-5)/p,
        sum(SP.est$act.coef.mat[,4]<sig_level)/p)
      temp_results_SIM[j, ] <- c(
        SIM.p.values[1],
        sum(abs(SIM.est$phi_est[-1])>1e-5)/p,
        sum(SIM.p.values[-1]<sig_level)/p)
    }
    results[i, ] <- colMeans(temp_results)
    results_SP[i, ] <- colMeans(temp_results_SP)
    results_SIM[i, ] <- colMeans(temp_results_SIM)
    
  }
}

# output perturbation results
disturb_df <- data.frame(
  degree = disturb_degree,
  p_value = results[, 1],
  active = results[, 2],
  significant = results[, 3]
)
write.csv(disturb_df, paste0(prefix, "disturb.csv"), row.names = FALSE)

disturb_df <- data.frame(
  degree = disturb_degree,
  p_value = results_SP[, 1],
  active = results_SP[, 2],
  significant = results_SP[, 3]
)
write.csv(disturb_df, paste0(prefix, "disturb SP.csv"), row.names = FALSE)

disturb_df <- data.frame(
  degree = disturb_degree,
  p_value = results_SIM[, 1],
  active = results_SIM[, 2],
  significant = results_SIM[, 3]
)
write.csv(disturb_df, paste0(prefix, "disturb SIM.csv"), row.names = FALSE)


# descriptive statistics
max(rowSums(A))
mean(rowSums(A))
max(colSums(A))
mean(colSums(A))


