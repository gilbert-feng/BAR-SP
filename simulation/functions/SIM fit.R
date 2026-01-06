# function to implement SIM estimates
SIMCV <- function(X, y, P, homo) {
  # X: covariate matrix
  # Y: response variable
  # P: the observed linking probability matrix
  # homo: 1: homoskedastic disturbance 0: heteroskedastic disturbance
  
  n <- nrow(X)
  p <- ncol(X)
  weight_total <- list(P)
  
  # initial estimates
  init <- est_initial(n, p, X, y, 1, weight_total)
  delta <- init$phi
  Wy <- init$Wy
  
  if (homo == 1) {
    # QMLE only use initial spatial parameters
    result <- qmle(y, Wy, X, delta[1], weight_total)
  } else {
    # GMM use all initial parameters
    result <- gmm(y, Wy, X, delta, weight_total)
  }
  
  phi_est <- result$phi_recover
  phi_sd <- result$se_recover
  
  # calculate mean squared error (MSE)
  total_mse <- tryCatch({
    Wn <- phi_est[1] * P
    X_coef <- phi_est[-1]
    y_pred <- solve(Diagonal(n) - Wn, X %*% X_coef)
    mean((y - y_pred)^2)
  }, error = function(e) {
    warning("Fail to calculate mean squared error: ", e$message)
    NA
  })
  
  list(
    phi_est = phi_est,
    phi_sd = phi_sd,
    in_mse = total_mse
  )
}

# calculate trace of A*B
tr_AB <- function(A, B) {
  if (ncol(A) != nrow(B)) {
    stop("Columns of A should equal to rows of B")
  }
  sum(A * t(B)) 
}

# initial estimate
est_initial <- function(n, p, X, y, weight_mat_num, weight_total) {
  # initialize Wy and WX
  Wy <- matrix(0, n, weight_mat_num)
  WX <- matrix(0, n, p * weight_mat_num * 2)
  
  # calculate weighted matrix
  for (i in 1:weight_mat_num) {
    start_idx <- (i - 1) * p + 1
    end_idx <- i * p
    
    # calculate WX
    WX[, start_idx:end_idx] <- weight_total[[i]] %*% X
    WX[, p * weight_mat_num + start_idx:end_idx] <- weight_total[[i]] %*% WX[, start_idx:end_idx]
    
    # calculate Wy
    Wy[, i] <- as.vector(weight_total[[i]] %*% y)  
  }
  
  # construct IV matrix
  Q <- cbind(X, WX)
  
  QZ <- t(Q) %*% cbind(Wy, X)
  QQ <- t(Q) %*% Q
  Qy <- t(Q) %*% y
  
  # make all components be matrix
  if (!is.matrix(QZ)) QZ <- as.matrix(QZ)
  if (!is.matrix(QQ)) QQ <- as.matrix(QQ)
  if (!is.matrix(Qy)) Qy <- as.matrix(Qy)
  
  phi <- solve(t(QZ) %*% solve(QQ, QZ), t(QZ) %*% solve(QQ, Qy))
  phi <- as.vector(phi)  
  
  list(phi = phi, Wy = Wy)
}

# GMM estimator 
gmm <- function(y, Wy, X, para, weight_total) {
  weight_mat_num <- length(weight_total)
  n <- nrow(X)
  p <- ncol(X)
  param_num <- p + weight_mat_num
  
  # define GMM objective function
  GMM <- function(para) {
    Sy <- y
    Wn <- Matrix(0, n, n, sparse = TRUE)
    
    if (length(para) != param_num) {
      stop("Incorrect parameter length")
    }
    
    for (i in 1:weight_mat_num) {
      Sy <- Sy - para[i] * as.vector(Wy[, i]) 
      Wn <- Wn + para[i] * weight_total[[i]]
    }
    
    S <- Diagonal(n) - Wn
    V <- Sy - X %*% para[(weight_mat_num + 1):length(para)]
    
    # initialize f
    f <- numeric(weight_mat_num + p)
    
    for (i in 1:weight_mat_num) {
      w_i_s <- weight_total[[i]] %*% solve(S)
      diag_term <- diag(diag(w_i_s))
      
      term1 <- as.numeric(t(V) %*% weight_total[[i]] %*% y)
      term2 <- as.numeric(t(V) %*% diag_term %*% V)
      f[i] <- term1 - term2
    }
    
    f[(weight_mat_num + 1):(weight_mat_num + p)] <- as.vector(t(X) %*% V)
    f
  }
  
  sol <- tryCatch(
    nleqslv(para, GMM, control = list(btol = 1e-3)),
    error = function(e) {
      warning("fail to solve the equation: ", e$message)
      list(x = rep(NA, param_num))
    }
  )
  
  if (any(is.na(sol$x))) {
    return(list(
      phi_recover = rep(NA, param_num),
      se_recover = rep(NA, param_num),
      Wn = Matrix(0, n, n),
      SE = matrix(NA, param_num, param_num)
    ))
  }
  
  delta <- sol$x
  
  # recover parameters
  lambda_hat <- delta[1:weight_mat_num]
  beta <- delta[(weight_mat_num + 1):length(delta)]
  phi_hat <- c(lambda_hat, beta)
  
  # calculate standard error
  Sy <- y
  Wn <- Matrix(0, n, n, sparse = TRUE)
  
  for (i in 1:weight_mat_num) {
    Sy <- Sy - lambda_hat[i] * as.vector(Wy[, i]) 
    Wn <- Wn + lambda_hat[i] * weight_total[[i]]
  }
  
  s <- Diagonal(n) - Wn
  Xb <- X %*% beta
  v_hat <- Sy - Xb
  Sigma_hat <- Diagonal(x = as.vector(v_hat)^2)
  
  # initialize Gamma and Omega
  Gamma <- matrix(0, param_num, param_num)
  Omega <- matrix(0, param_num, param_num)
  
  for (i in 1:weight_mat_num) {
    w_i_s <- weight_total[[i]] %*% solve(s)
    w_i_sXb <- as.vector(w_i_s %*% Xb)  
    Tg_i <- w_i_s - diag(diag(w_i_s))
    
    for (j in i:weight_mat_num) {
      w_j_s <- weight_total[[j]] %*% solve(s)
      w_j_sXb <- as.vector(w_j_s %*% Xb)  
      Tg_j <- w_j_s - diag(diag(w_j_s))
      
      Gamma[j, i] <- tr_AB(t(w_i_s) %*% (Tg_j + t(Tg_j)), as.matrix(Sigma_hat)) + 
        sum(w_i_sXb * w_j_sXb)  
      
      Omega[j, i] <- tr_AB(as.matrix(Tg_j %*% Sigma_hat), as.matrix((Tg_i + t(Tg_i)) %*% Sigma_hat)) +
        sum(w_j_sXb * as.vector(Sigma_hat %*% w_i_sXb))  
    }
    
    # Gamma matrix's beta component 
    Gamma[(weight_mat_num + 1):param_num, i] <- as.vector(t(X) %*% w_i_sXb)
    
    # Omega matrix's beta component
    Omega[(weight_mat_num + 1):param_num, i] <- as.vector(t(X) %*% (Sigma_hat %*% w_i_sXb))
  }
  
  # fill up lower triangle elements
  Gamma <- Gamma + t(Gamma) - diag(diag(Gamma))
  Omega <- Omega + t(Omega) - diag(diag(Omega))
  
  # beta component
  Omega[(weight_mat_num + 1):param_num, (weight_mat_num + 1):param_num] <- as.matrix(t(X) %*% Sigma_hat %*% X)
  Gamma[(weight_mat_num + 1):param_num, (weight_mat_num + 1):param_num] <- as.matrix(t(X) %*% X)
  
  # calculate covariance matrix
  SE <- tryCatch(
    solve(Gamma) %*% Omega %*% solve(t(Gamma)),
    error = function(e) {
      warning("Fail to calculate covariance matrix: ", e$message)
      matrix(NA, param_num, param_num)
    }
  )
  
  list(
    phi_recover = phi_hat,
    se_recover = sqrt(abs(diag(SE))),
    Wn = Wn,
    SE = SE
  )
}

# QMLE estimator
qmle <- function(y, Wy, X, para, weight_total) {
  weight_mat_num <- length(weight_total)
  n <- nrow(X)
  p <- ncol(X)
  param_num <- p + weight_mat_num
  
  # define the concentrated likelihood function
  CMLE <- function(para) {
    Sy <- y
    Wn <- Matrix(0, n, n, sparse = TRUE)
    
    if (length(para) != weight_mat_num) {
      stop("Incorrect parameter length")
    }
    
    for (i in 1:weight_mat_num) {
      Sy <- Sy - para[i] * as.vector(Wy[, i])  
      Wn <- Wn + para[i] * weight_total[[i]]
    }
    
    beta <- solve(t(X) %*% X, t(X) %*% Sy)
    s <- Diagonal(n) - Wn
    
    log_det <- determinant(s, logarithm = TRUE)$modulus[1]
    
    v_hat <- Sy - X %*% beta
    sige <- mean(v_hat^2)
    
    n/2 * log(sige) - log_det
  }
  
  # optimize likelihood function
  opt <- tryCatch(
    optim(para, CMLE, method = "BFGS", control = list(maxit = 1000)),
    error = function(e) {
      warning("Fail to optimize: ", e$message)
      list(par = rep(NA, weight_mat_num), value = NA)
    }
  )
  
  if (any(is.na(opt$par))) {
    return(list(
      phi_recover = rep(NA, param_num),
      se_recover = rep(NA, param_num),
      Wn = Matrix(0, n, n),
      SE = matrix(NA, param_num, param_num)
    ))
  }
  
  lambda_hat <- opt$par
  
  # recover parameters
  Sy <- y
  for (i in 1:weight_mat_num) {
    Sy <- Sy - lambda_hat[i] * as.vector(Wy[, i])  
  }
  beta <- as.vector(solve(t(X) %*% X, t(X) %*% Sy))  
  phi_hat <- c(lambda_hat, beta)
  
  # calculate standard error
  Wn <- Matrix(0, n, n, sparse = TRUE)
  Sy <- y
  for (i in 1:weight_mat_num) {
    Sy <- Sy - lambda_hat[i] * as.vector(Wy[, i])  
    Wn <- Wn + lambda_hat[i] * weight_total[[i]]
  }
  
  s <- Diagonal(n) - Wn
  Xb <- X %*% beta
  v_hat <- Sy - Xb
  sige <- mean(v_hat^2)
  mu3 <- mean(v_hat^3)
  mu4 <- mean(v_hat^4)
  dev <- mu4 - 3 * sige^2
  
  # initialize Hessian and Omega
  Hessian <- matrix(0, param_num + 1, param_num + 1)
  Omega <- matrix(0, param_num + 1, param_num + 1)
  
  for (i in 1:weight_mat_num) {
    G_i <- weight_total[[i]] %*% solve(s)
    G_iXb <- as.vector(G_i %*% Xb)  
    
    for (j in i:weight_mat_num) {
      G_j <- weight_total[[j]] %*% solve(s)
      G_jXb <- as.vector(G_j %*% Xb)  
      
      term1 <- sum(G_iXb * G_jXb)  
      term2 <- sige * tr_AB(G_i + t(G_i), as.matrix(G_j))
      Hessian[j, i] <- term1 + term2
      
      if (j == i) {
        Omega[j, i] <- dev * sum(diag(G_i)^2) + 
          2 * mu3 * sum(G_iXb * diag(G_i))
      } else {
        Omega[j, i] <- dev * sum(diag(G_i) * diag(G_j)) + 
          mu3 * (sum(G_iXb * diag(G_j)) + 
                   sum(G_jXb * diag(G_i)))
      }
    }
    
    # Hessian matrix's beta component
    Hessian[(weight_mat_num + 1):param_num, i] <- as.vector(t(X) %*% G_iXb)
    
    # Hessian matrix's sigma component
    Hessian[param_num + 1, i] <- sum(diag(G_i))
    
    # Omega matrix's beta component
    Omega[(weight_mat_num + 1):param_num, i] <- mu3 * as.vector(t(X) %*% diag(G_i))
    
    # Omega matrix's sigma component
    Omega[param_num + 1, i] <- 0.5 * (dev * sum(diag(G_i)) + 
                                        mu3 * sum(G_iXb))
  }
  
  
  # beta and sigma components
  Hessian[(weight_mat_num + 1):param_num, (weight_mat_num + 1):param_num] <- as.matrix(t(X) %*% X)
  Hessian[param_num + 1, param_num + 1] <- n/(2 * sige)
  diag_Hessian <- diag(diag(Hessian))
  temp_Hessian <- Hessian - diag_Hessian
  Hessian <- Hessian + t(temp_Hessian)
  Hessian <- Hessian / sige  
  
  # fill up Omega matrix's remaining elements 
  X_colsums <- colSums(X) 
  Omega[param_num + 1, (weight_mat_num + 1):param_num] <- (mu3 / (2 * sige)) * X_colsums
  Omega[param_num + 1, param_num + 1] <- n * dev / (4 * sige^2)
  
  # symmetrization
  diag_Omega <- diag(diag(Omega))
  temp_Omega <- Omega - diag_Omega
  Omega <- Omega + t(temp_Omega)
  Omega <- Omega / (sige^2) 
   
   # calculate the covariance matrix
    SE <- tryCatch({
      H_inv <- solve(Hessian)
      H_inv %*% (Hessian + Omega) %*% H_inv
    }, error = function(e) {
      warning("Fail to calculate the covariance matrix: ", e$message)
      matrix(NA, param_num + 1, param_num + 1)
    })
   
   # remove sigma component 
   SE <- SE[1:param_num, 1:param_num]
   

   list(
     phi_recover = phi_hat,
     se_recover = sqrt(abs(diag(SE))),
     Wn = Wn,
     SE = SE
   )
}

