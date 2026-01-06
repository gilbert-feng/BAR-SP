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
n_seq <- 1000
n_degree <- 3
rep_time <- 100
G <- 4
K <- 4
r <- 1
sigma_0 <- 1
p_seq <- 1:10*10*2
s_0 <- 5

# True coefficients
gamma_0 <- rep(1,4)

# Hyper-parameters
iter_max <- 1e3
conv_crit <- 1e-5
sig_level <- 0.05
K_fold <- 2
lambda_seq <- NULL
#SBM_type <- "sample" # sample/MLE
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
temp_mse_mat <- matrix(0, rep_time, 2)
mse_mat <- matrix(0,length(p_seq),3)

# Main loop
for (homo in c(0, 1)) {
  for (SBM_type in c("sample", "MLE")) {
    for (n in n_seq) {
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
        eigen.P <- eigs_sym(A=P, k=4)
        U <- eigen.P$vectors[,1:4]
        alpha.coef <- matrix(sqrt(n) * gamma_0, ncol=1)
        alpha <- U %*% alpha.coef
        
        for (p_ind in 1:length(p_seq)) {
          p <- p_seq[p_ind]
          # X generation
          X <- matrix(rnorm(n*p),n,p)
          G_num <- p/G
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
          beta_0 <- rep(0,p)
          beta_0[2:(2+s_0)] <- 1
          beta_0 <- matrix(beta_0, ncol=1)
          Xbeta <- X %*% beta_0
          EY <- Xbeta + alpha
          
          # Calculate true projected components
          ProdSVD.fit <- ProdSVD(X,P,K,r=r) 
          SP.true <- SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,EY,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,homo,sig_level)
          SP.true.beta <- SP.true$beta
          
          temp_X_proj <- ProdSVD.fit$x_proj
          temp_Z <- ProdSVD.fit$Z
          
          
          cat("------ Simulation start for case: homo -", homo, "n =", n, "\n")
          
          if (homo == 1) {
            prefix <- paste0("results G"," mse p","/homo/")
          } else {
            prefix <- paste0("results G"," mse p","/hete/")
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
            # Estimate product SVD
            A <- net.gen.from.P(P)
            if (SBM_type == "MLE") {
              SBM.MLE <- SBM.estimate(A,g=big.model$g)
              A <- SBM.MLE$Phat
            }
            eigen_A <- eigs_sym(A=A,k=K)
            W <- matrix(eigen_A$vectors[,1:K],ncol=K)
            temp_ProdSVD <- svd(t(temp_Z)%*%W,nu=p,nv=K)
            S_hat <- temp_ProdSVD$d
            Z_hat <- temp_Z%*%temp_ProdSVD$u
            W_cup <- W%*%temp_ProdSVD$v
            ProdSVD.fit <- list(x_proj = temp_X_proj, S_hat = S_hat,
                                Z = temp_Z, Z_hat = Z_hat,
                                W = W, W_cup = W_cup,
                                r = r)
            
            
            # (1) A-SP method
            start_time <- Sys.time()
            BAR.SP.est <- BAR.SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,homo,K_fold,lambda_seq,sig_level)
            temp_mse_mat[i, 1] <- sqrt(mean((BAR.SP.est$beta[-1] - SP.true.beta[-1])^2))
            cat("Simulation process (adaptive)", round(100 * i / rep_time), "% complete!\n")
            
            # (2) SP 
            start_time <- Sys.time()
            SP.est <- SP.Prod(X,ProdSVD.fit$x_proj,ProdSVD.fit$Z_hat,Y,ProdSVD.fit$W_cup,K,ProdSVD.fit$r,homo,sig_level)
            temp_mse_mat[i, 2] <- sqrt(mean((SP.est$beta[-1] - SP.true.beta[-1])^2))
            cat("Simulation process (SP)", round(100 * i / rep_time), "% complete!\n")
            
            
          }
          mse_mat[p_ind,degree] <- mean(temp_mse_mat[,1]/temp_mse_mat[,2])
          
          # simulation dialogue
          cat("------ Simulation end for case: homo -", homo, "n =", n, "\n")
          
          
        }
        
        
        }
    }
    write.csv(mse_mat, paste0(prefix, SBM_type, " mse.csv"))
    
  }
  
}