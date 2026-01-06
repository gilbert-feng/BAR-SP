# function to implement regression with network cohesion (RNC) method
RNCCV <- function(X, y, P, lambda_seq = NA, gamma, K_fold) {
  # X: covariate matrix
  # Y: response variable
  # P: the observed linking probability matrix
  # lambda_seq: the hyper-parameter sequence
  # gamma: a small positive constant to obtain a regularized Laplacian
  # K_fold: fold number for cross-validation
  
  n <- nrow(X)
  p <- ncol(X)
  Xy_max <- max(abs(t(X) %*% y))
  
  if (is.na(lambda_seq[1])) {
    # lambda_seq <- logspace(log10(1e-3), log10(Xy_max), 10)
    lambda_seq <- seq(1, Xy_max, length.out = 10) * 10
  }
  
  L <- diag(rowSums(P)) - P
  
  XX <- t(X) %*% X
  X_tild <- cbind(diag(n), X)
  
  pred_mse_mat <- matrix(0, length(lambda_seq), K_fold)
  
  # sample split and hyper-parameter selection
  groups <- sample(rep(1:K_fold, length.out = n))
  for (fold_id in 1:K_fold) {
    test_index <- groups == fold_id
    train_index <- !test_index
    X_train <- X[train_index, ]
    n_tr <- nrow(X_train)
    X_test <- X[test_index, ]
    XX_train <- t(X_train) %*% X_train
    X_train_tild <- cbind(diag(n_tr), X_train)
    
    for (i in 1:length(lambda_seq)) {
      L_train <- L[train_index, train_index]
      XX_train_tilde <- rbind(
        cbind(diag(n_tr) + lambda_seq[i] * (L_train + gamma * diag(n_tr)), X_train),
        cbind(t(X_train), XX_train)
      )
      theta <- solve(XX_train_tilde) %*% t(X_train_tild) %*% y[train_index]
      beta_train <- theta[(length(theta) - p + 1):length(theta)]
      alpha_train <- theta[1:(length(theta) - p)]
      
      L_test <- L[test_index, test_index]
      L_21 <- L[test_index, train_index]
      
      # add small perturbation to diagnal elements 
      alpha_test <- -solve(L_test + 1e-5 * diag(nrow(L_test))) %*% L_21 %*% alpha_train
      pred_mse_mat[i, fold_id] <- mean((y[test_index] - X_test %*% beta_train - alpha_test)^2)
    }
  }
  
  pred_mse <- rowMeans(pred_mse_mat)
  if (length(lambda_seq) > 1) {
    min_index <- which.min(pred_mse)
    min_mse <- pred_mse[min_index]
    lambda_opt <- lambda_seq[min_index]
  } else {
    min_mse <- min(pred_mse)
    lambda_opt <- lambda_seq
  }
  
  XX_tilde <- rbind(
    cbind(diag(n) + lambda_opt * (L + gamma * diag(n)), X),
    cbind(t(X), XX)
  )
  theta <- solve(XX_tilde) %*% t(X_tild) %*% y
  
  beta_hat <- theta[(length(theta) - p + 1):length(theta)]
  alpha_hat <- theta[1:(length(theta) - p)]
  
  total_mse <- mean((y - (X %*% beta_hat + alpha_hat))^2)
  
  list(beta_hat = beta_hat, alpha_hat = alpha_hat, min_mse = min_mse, 
       in_mse = total_mse, lambda_opt = lambda_opt)
}    