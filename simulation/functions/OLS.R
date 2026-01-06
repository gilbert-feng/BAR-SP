# function to implement OLS estimation
OLSCV <- function(X, y, homo) {
  # X: covariate matrix
  # Y: response variable
  # homo: 1: homoskedastic disturbance 0: heteroskedastic disturbance
  
  n <- nrow(X)
  p <- ncol(X)
  XX <- t(X) %*% X
  # XX_inv = solve(XX)
  beta_hat <- solve(XX) %*% t(X) %*% y
  
  # variance estimate
  if (homo == 1) {
    # sigma2_hat = sum((y - (PC_y +  xi_hat + alpha_hat))^2)/(n-p_X-K+r); % for homo, this is better
    sigma2_hat <- sum((y - (X %*% beta_hat))^2)/(n - p)
    beta_sd <- sigma2_hat * solve(XX)
  } else {
    # D_e = diag(y - (PC_y +  xi_hat + alpha_hat));
    # D_e = diag(y - (X*beta_hat +  xi_hat + alpha_hat)); % better
    # D_e = diag(y - (X*(beta_hat +  theta_hat) + alpha_hat)); % similar to
    # better
    D_e <- diag(as.vector(y - (X %*% beta_hat)))
    Sigma <- solve(XX) %*% t(X) %*% D_e
    beta_sd <- Sigma %*% t(Sigma)
  }
  
  total_mse <- mean((y - (X %*% beta_hat))^2)
  
  list(beta_hat = beta_hat, beta_sd = beta_sd, in_mse = total_mse)
}
