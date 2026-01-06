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

# Seed
set.seed(1234)

# Parameters
n_seq <- c(3, 5, 10, 20, 40) * 100
n_degree <- 3
rep_time <- 800
prob_in <- 1
prob_out <- 0.2
G <- 4
K <- 4
r <- 1
p_0 <- 12
p <- p_0 + r
sigma_0 <- 1

# True coefficients
beta_0 <- rep(0,p)
beta_0[2:7] <- 1
beta_0 <- matrix(beta_0, ncol=1)
#beta_0 <- matrix(c(0, rep(1, 6), rep(0, p_0 - 6)), ncol=1)

# Hyper-parameters
iter_max <- 1e3
conv_crit <- 1e-5
sig_level <- 0.05
K_fold <- 2
lambda_seq <- NULL
SBM_type <- "MLE" # sample/MLE
# row_names <- c('bias(SP)', 'rmse(SP)', 'std(SP)', 'se(SP)', 'diff(SP)',
#                'bias(ridge)', 'rmse(ridge)', 'std(ridge)', 'se(ridge)', 'diff(ridge)',
#                'bias(adaptive)', 'rmse(adaptive)', 'std(adaptive)', 'se(adaptive)', 'diff(adaptive)')

# Covariance structure
# rho <- 0.95
# C <- matrix(0, p, p)
# for (i in 1:p) {
#   for (j in 1:p) {
#     C[i, j] <- rho^abs(i - j)
#   }
# }
# C <- chol(C)

# Storage initialization
total_param_num <- p
param_hat <- array(0, dim = c(rep_time, total_param_num, 2))
SP_beta_true <- matrix(0,p,5) # beta_0 SP_beta_0 SP BAR 
SP_beta_true[,1] <- beta_0
temp_lm_mat <- matrix(0,rep_time,p)
param_hat_sd <- array(0, dim = c(rep_time, total_param_num, 2))
measures_num <- 6 # project measures + network p-value
temp_measures <- matrix(0,rep_time,measures_num)
measures <- matrix(0, length(n_seq)*measures_num, 3) 
time_mat <- matrix(0, rep_time, 2)
zero_num <- 2
temp_zero <- array(0, dim = c(rep_time, 2, 3))
zero_mat <- matrix(0, length(n_seq)*zero_num, 9)
sd_ratio_table <- array(0, dim = c(length(n_seq), n_degree * 3, 2))

for (SBM_type in c("sample", "MLE")) {
  for (gamma_effect in c(0,1)) {
    gamma_0 <- rep(gamma_effect,4)
    # Main loop
    for (homo in c(0,1)) {
      n_count <- 0
      for (n in n_seq) {
        n_count <- n_count + 1
        for (degree in 1:n_degree) {
          # network generation
          if (degree == 1) {
            big.model <- BlockModel.Gen(lambda=2*log(n), n=n, beta=0.2, K=K)
          } else if (degree == 2) {
            big.model <- BlockModel.Gen(lambda=sqrt(n), n=n, beta=0.2, K=K)
          } else {
            big.model <- BlockModel.Gen(lambda=n^(2/3), n=n, beta=0.2, K=K)
          }
          P <- big.model$P
          A <- net.gen.from.P(P)
          if (SBM_type == "MLE") {
            SBM.MLE <- SBM.estimate(A,g=big.model$g)
            A <- SBM.MLE$Phat
          }
          
          eigen.P <- eigs_sym(A=P, k=4)
          U <- eigen.P$vectors[,1:4]
          alpha.coef <- matrix(sqrt(n) * gamma_0, ncol=1)
          alpha <- U %*% alpha.coef
          
          # X generation
          X <- matrix(rnorm(n*p),n,p)
          G_num <- p_0/G
          for (i in 1:G) {
            start_index <- (i-1)*G_num + 1
            X[,start_index+r] <- rnorm(n,0,1)
            for (j in start_index + 1:2) {
              X[,j+r] <- X[,start_index+r]+(runif(n)-0.5)*4
            }
          }
          
          if (r>=1) {
            X[,1:r] <- U[,1:r]
          }
          X <- scale(X, center=TRUE, scale=TRUE) * sqrt(n/(n-1))
          Xbeta <- X %*% beta_0
          EY <- Xbeta + alpha
          
          
          # Calculate true projected components
          ProdSVD.fit <- ProdSVD(X,P,K,r=r) 
          SP.true <- SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,EY,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,homo,sig_level)
          SP.true.beta <- SP.true$beta
          SP_beta_true[,2] <- SP.true.beta
          SP_beta_active_ind <- which(beta_0!=0)
          #SP_beta_active_ind <- which(abs(SP.true.beta)>1e-5)
          
          
          # Estimate product SVD
          ProdSVD.fit <- ProdSVD(X,A,K,r=r)
          
          
          cat("------ Simulation start for case: homo -", homo, "n =", n, "\n")
          
          if (homo == 1) {
            prefix <- paste0("results G ", "gamma = ",gamma_effect, " SBM_type = ", SBM_type, "/homo/")
          } else {
            prefix <- paste0("results G ", "gamma = ",gamma_effect, " SBM_type = ", SBM_type, "/hete/")
          }
          dir.create(paste0(prefix, "n = ", n), recursive = TRUE, showWarnings = FALSE)
          
          for (i in 1:rep_time) {
            if (homo == 1) {
              epsilon <- rnorm(n, 0, sigma_0)
            } else {
              temp_sd <- sqrt(1/2+1/2*(X[, 2])^2)
              epsilon <- rnorm(n, 0, sigma_0) * temp_sd
            }
            
            Y <- EY + epsilon
            
            lm.fit <- lm(Y~X-1)
            temp_lm_mat[i,] <- lm.fit$coefficients
            # (1) A-SP method
            start_time <- Sys.time()
            BAR.SP.est <- BAR.SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,homo,K_fold,lambda_seq,sig_level)
            time_mat[i, 1] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
            active_ind <- BAR.SP.est$active_ind
            param_hat[i, , 1] <- BAR.SP.est$beta
            temp_cov_sd <- rep(0,p)
            temp_cov_sd[active_ind] <- sqrt(diag(BAR.SP.est$act_cov_hat))
            param_hat_sd[i, , 1] <- temp_cov_sd
            temp_measures[i,] <- c(BAR.SP.est$project_measures,BAR.SP.est$chisq.p)
            active_model <- ZeroCount(SP.true.beta[-1],BAR.SP.est$beta[-1],1e-5);
            G_index <- GroupCount(BAR.SP.est$beta, G, r, 1e-5)
            temp_zero[i,1,1] <- active_model$Precision
            temp_zero[i,1,2] <- active_model$Recall
            temp_zero[i,1,3] <- G_index
            
            cat("Simulation process (adaptive)", round(100 * i / rep_time), "% complete!\n")
            
            # (2) SP 
            start_time <- Sys.time()
            SP.est <- SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,homo,sig_level)
            time_mat[i, 2] <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
            param_hat[i, , 2] <- SP.est$beta
            param_hat_sd[i, , 2] <- sqrt(diag(SP.est$act_cov_hat))
            active_model <- ZeroCount(SP.true.beta[-1],SP.est$beta[-1],1e-5);
            G_index <- GroupCount(SP.est$beta, G, r, 1e-5)
            temp_zero[i,2,1] <- active_model$Precision
            temp_zero[i,2,2] <- active_model$Recall
            temp_zero[i,2,3] <- G_index
            
            cat("Simulation process (SP)", round(100 * i / rep_time), "% complete!\n")
            
            
          }
          # SP_beta_true[,3] <- colMeans(param_hat[, , 1])
          # SP_beta_true[,4] <- colMeans(param_hat[, , 2])
          # SP_beta_true[,5] <- colMeans(temp_lm_mat)
          SP_beta_true[,3] <- param_hat[1, , 1]
          SP_beta_true[,4] <- param_hat[1, , 2]
          SP_beta_true[,5] <- temp_lm_mat[1,]
          measures[((n_count-1)*measures_num+1):(n_count*measures_num),degree] <- colMeans(temp_measures)
          temp_zero_mean <- apply(temp_zero,2:3,mean)
          zero_mat[((n_count-1)*zero_num+1):(n_count*zero_num),degree] <- temp_zero_mean[,1]
          zero_mat[((n_count-1)*zero_num+1):(n_count*zero_num),degree+3] <- temp_zero_mean[,2]
          zero_mat[((n_count-1)*zero_num+1):(n_count*zero_num),degree+6] <- temp_zero_mean[,3]
          
          for (model_id in 1:2) {
            # rmse for all param
            param_rmse <- sqrt(colMeans((param_hat[, , model_id] - matrix(SP.true.beta, nrow = rep_time, ncol = length(SP.true.beta), byrow = TRUE))^2))
            # CI for oracle active
            conf_var <- 1.96 * param_hat_sd[, SP_beta_active_ind, model_id]  # 注意索引偏移
            conf_u <- param_hat[, SP_beta_active_ind, model_id] + conf_var
            conf_l <- param_hat[, SP_beta_active_ind, model_id] - conf_var
            conv_prob <- mean((t(conf_u) > SP.true.beta[SP_beta_active_ind]) & (t(conf_l) < SP.true.beta[SP_beta_active_ind]))
            conf_length <- mean(2 * conf_var)
            bs_rmse <- mean(param_rmse[SP_beta_active_ind])
            
            sd_ratio_table[n_count, degree + c(0, n_degree, 2 * n_degree), model_id] <- 
              c(conf_length, conv_prob, bs_rmse)
          }
          
          # Save results
          write.csv(SP_beta_true, paste0(prefix, "n = ", n, "/SP_ave_beta degree_", degree, ".csv"))
          
          # simulation diaglogue
          cat("------ Simulation end for case: homo -", homo, "n =", n, "\n")
        }
      }
      
      write.csv(measures, paste0(prefix, "measures.csv"))
      write.csv(zero_mat, paste0(prefix, "zeros.csv"))
      
      
      # relative measure
      for (model_id in 1) {
        # Relative RMSE
        sd_ratio_table[, (n_degree * 2 + 1):(n_degree * 3), model_id] <- 
          sd_ratio_table[, (n_degree * 2 + 1):(n_degree * 3), model_id] / 
          sd_ratio_table[, (n_degree * 2 + 1):(n_degree * 3), 2]
        
        # Relative CI length
        sd_ratio_table[, 1:n_degree, model_id] <- 
          sd_ratio_table[, 1:n_degree, model_id] / sd_ratio_table[, 1:n_degree, 2]
      }
      
      # save results
      for (model_id in 1:2) {
        write.csv(sd_ratio_table[, , model_id], paste0(prefix, "model_", model_id, "_sd_ratio.csv"))
      }
    }
  }
  
}

