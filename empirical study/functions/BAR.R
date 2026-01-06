
# preliminary function to implement BAR estimation 
BAR.SP.CV <- function(X,Y,A,K,homo,lambda_seq,K_fold=2,r=NULL,thr=NULL,alpha.CI=0.05,boot.thr=TRUE,boot.n=50, regular = 0, cv_sqrt = 0) {
  # X: covariate matrix
  # Y: response variable
  # A: the observed linking probability matrix
  # K: the dimension of network spectral information
  # homo: 1: homoskedastic disturbance 0: heteroskedastic disturbance
  # lambda_seq: BAR regularization hyper-parameter
  # K_fold: fold number for cross-validation
  # r: dimension of common component
  # thr: threshold for determing the dimension of common component
  # alpha.CI: significance level
  # boot.thr: whether to use bootstrap threshold to determine the dimension of common component
  # boot.n: the iteration number of bootstrap for for determing the dimension of common component
  # regular: whether add small perturbation to make OLS stable
  # cv_sqrt: whether use sqrt_n regularization
  
  # Subpace projection
  ProdSVD <- ProdSVD(X,A,K,r=r,thr=thr,boot.thr=boot.thr,boot.n=boot.n)
  x_proj <- ProdSVD$x_proj
  Z_hat <- ProdSVD$Z_hat
  W_cup <- ProdSVD$W_cup
  r <- ProdSVD$r
  thr <- ProdSVD$thr
  # BAR estimation
  BAR.SP.est <- BAR.SP.Prod(X,x_proj,Z_hat,Y,W_cup,K,r,homo,K_fold,lambda_seq,alpha.CI, regular = regular, cv_sqrt = cv_sqrt)

  return(c(BAR.SP.est,r=r,thr=thr))
}



# intermediate  function to implement BAR estimation (free of repeated SVD)
BAR.SP.Prod <- function(X,x_proj,Z_hat,Y,W_cup,K,r,homo,K_fold=2,lambda_seq,alpha.CI,active_thr=1e-5, regular = 0, cv_sqrt = 0) {
  # X: covariate matrix
  # x_proj: covariate-related projection matrix
  # Z_hat: rotated orthonormal covariate matrix
  # lambda_seq: BAR regularization hyper-parameter
  # Y: response variable
  # W_cup: rotated network spectral matrix
  # K: the dimension of network spectral information
  # r: dimension of common component
  # homo: 1: homoskedastic disturbance 0: heteroskedastic disturbance
  # K_fold: fold number for cross-validation
  # lambda_seq: BAR regularization hyper-parameter
  # alpha.CI: significance level
  # active_thr: the threshold for determing whether a coefficient is 0 or not
  # regular: whether add small perturbation to make OLS stable
  # cv_sqrt: whether use sqrt_n regularization
  
  if (is.null(lambda_seq)) {
    Xy_max <- max(abs(t(X)%*%Y))
    lambda_seq <- 10^seq(log10(1e-3), log10(Xy_max), length.out = 10)
  }
  # I Projection
  if(r > 0){
    M <- cbind(Z_hat[,-(1:r)],W_cup[,-(1:r)])
    P_CN <- M%*%solve(t(M)%*%M,t(M))
    P_R <- matrix(Z_hat[,1:r],ncol=r)%*%t(matrix(Z_hat[,1:r],ncol=r))
    P_C <- cbind(Z_hat[,-(1:r)],matrix(0,nrow=n,ncol=K-r))%*%solve(t(M)%*%M,t(M))
    P_N <- P_CN - P_C
    H <- P_CN + P_R
    theta_hat <- x_proj%*%P_R%*%Y
    xi_hat <- P_R%*%Y
    PU.hat <- t(W_cup[,-(1:r)])
  }else{
    M <- cbind(Z_hat,W_cup)
    P_CN <- M%*%solve(t(M)%*%M,t(M))
    xi_hat <- matrix(0,nrow(X),1)
    P_C <- cbind(Z_hat,matrix(0,nrow=n,ncol=K))%*%solve(t(M)%*%M,t(M))
    P_N <- P_CN - P_C
    H <- P_CN
    theta_hat <- matrix(0,nrow=ncol(X),ncol=1)
    PU.hat <- t(W_cup)
  }
  # II Recover parameters and measures
  beta_z <- P_C%*%Y
  alpha_hat <- P_N%*%Y
  fitted.Y <- H%*%Y
  project_measures <- matrix(rep(0,5),1)
  colnames(project_measures) <- c("P_C", "P_N", "P_CN", "P_R", "Error")
  project_measures[1] <- sum((beta_z)^2)
  project_measures[2] <- sum((alpha_hat)^2)
  project_measures[3] <- 2*sum(beta_z*alpha_hat)
  project_measures[4] <- sum((xi_hat)^2)
  project_measures <- project_measures/sum((Y)^2)
  project_measures[5] <- 1 - sum(project_measures[1:4])
  
  
  # III BAR estimation
  # CV
  if (length(lambda_seq) > 1) {
    pred_mse_mat <- matrix(0, nrow = length(lambda_seq), ncol = K_fold)
    #out_mse_mat <- matrix(0, nrow = length(lambda_seq), ncol = K_fold)
    # Generate K-fold groups
    groups <- sample(rep(1:K_fold, length.out = n))
    
    for (fold_id in 1:K_fold) {
      test_index <- (groups == fold_id)
      train_index <- !test_index
      X_train <- X[train_index, ]
      X_test <- X[test_index, ]
      
      # Training pseudo response
      PC_y_train <- P_C[train_index, train_index] %*% Y[train_index]
      
      for (i in seq_along(lambda_seq)) {
        # BAR iteration estimator
        if (cv_sqrt == 0) {
          beta_new <- BAR.iter(X_train,PC_y_train,lambda_seq[i],K_fold, regular = regular)
          
        } else {
          beta_new <- BAR.iter(X_train,PC_y_train,sqrt(nrow(X_train))*lambda_seq[i],K_fold, regular = regular)
          
        }
        # Prediction MSE
        y_pred <- X_test %*% beta_new
        y_proj <- P_C[test_index, test_index] %*% Y[test_index]
        pred_mse_mat[i, fold_id] <- mean((y_proj - y_pred)^2)
        #out_mse_mat[i, fold_id] <- mean((Y[test_index]-y_pred-xi_hat[test_index]-alpha_hat[test_index])^2)
      }
    }
    
    # Choose optimal lambda
    pred_mse <- rowMeans(pred_mse_mat)
    #out_mse <- rowMeans(out_mse_mat)
    min_index <- which.min(pred_mse)
    lambda_opt <- lambda_seq[min_index]
    
    min_mse <- pred_mse[min_index]
    #out_min_mse <- out_mse[min_index]
    
    
  } else {
    # Choose specified lambda
    lambda_opt <- lambda_seq
    min_mse <- NULL
  }
  
  # BAR estimate with optimal lambda
  beta_hat <- BAR.iter(X,beta_z,lambda_opt,K_fold, regular = regular)
  beta_hat[abs(beta_hat)<active_thr] <- 0
  active_ind <- which(abs(beta_hat)>active_thr)
  beta_act <- beta_hat[active_ind]
  
  # IV Inference
  ## (1) homo
  if (homo == 1) {
    sigma2_hat <- sum((Y-H%*%Y)^2)/(n-ncol(X)-K+r)
    
    # network inference
    cov.gamma <- PU.hat%*%P_N
    cov.gamma <- sigma2_hat*tcrossprod(cov.gamma)
    gamma <- PU.hat%*%alpha_hat
    if((K-r)>1){
      eig.cov.gamma <- eigen(cov.gamma,symmetric=TRUE)
      inv.sqrt.cov.gamma <- t(eig.cov.gamma$vectors)*sqrt(1/eig.cov.gamma$values)
      adjusted.gamma <- inv.sqrt.cov.gamma%*%gamma
    }else{
      adjusted.gamma <- 0
    }
    chisq.val <- sum(adjusted.gamma^2)
    chisq.p <- pchisq(chisq.val,df=K-r,lower.tail=FALSE)
    
    # BAR inference with active components
    tmp_cov <- t(X)%*%P_C
    tmp_cov <- tcrossprod(tmp_cov)
    X_active <- X[,active_ind]
    XX_active <- crossprod(X_active)
    cov_hat <- sigma2_hat*solve(XX_active, tmp_cov[active_ind, active_ind]) %*% solve(XX_active)
    diag_sd <- sqrt(diag(cov_hat))
    abs.t.val <- abs(beta_act/diag_sd)
    t.pval <- 1-pnorm(abs.t.val)
    CI.lower <- beta_act - qnorm(1-alpha.CI/2)*diag_sd
    CI.upper <- beta_act + qnorm(1-alpha.CI/2)*diag_sd
    coef.mat <- cbind(beta_act,CI.lower,CI.upper,t.pval)
    colnames(coef.mat) <- c("active coef","CI-lower","CI-upper","p-val")
  } else {
    sigma2_hat = NULL
    D_e <- Y - alpha_hat - xi_hat - X%*%beta_hat
    
    # network inference
    cov.gamma <- PU.hat%*%P_N%*%diag(as.vector(D_e))
    cov.gamma <- tcrossprod(cov.gamma)
    gamma <- PU.hat%*%alpha_hat
    if((K-r)>1){
      eig.cov.gamma <- eigen(cov.gamma,symmetric=TRUE)
      inv.sqrt.cov.gamma <- t(eig.cov.gamma$vectors)*sqrt(1/eig.cov.gamma$values)
      adjusted.gamma <- inv.sqrt.cov.gamma%*%gamma
    }else{
      adjusted.gamma <- 0
    }
    chisq.val <- sum(adjusted.gamma^2)
    chisq.p <- pchisq(chisq.val,df=K-r,lower.tail=FALSE)
    
    # BAR inference with active components
    tmp_cov <- t(X)%*%P_C%*%diag(as.vector(D_e))
    tmp_cov <- tcrossprod(tmp_cov)
    X_active <- X[,active_ind]
    XX_active <- crossprod(X_active)
    cov_hat <- solve(XX_active, tmp_cov[active_ind, active_ind]) %*% solve(XX_active)
    diag_sd <- sqrt(diag(cov_hat))
    abs.t.val <- abs(beta_act/diag_sd)
    t.pval <- 1-pnorm(abs.t.val)
    CI.lower <- beta_act - qnorm(1-alpha.CI/2)*diag_sd
    CI.upper <- beta_act + qnorm(1-alpha.CI/2)*diag_sd
    coef.mat <- cbind(beta_act,CI.lower,CI.upper,t.pval)
    colnames(coef.mat) <- c("active coef","CI-lower","CI-upper","p-val")
    
  }
  
  rownames(coef.mat) <- paste("V",active_ind,sep="")
  
  in_mse <- mean((Y-X%*%beta_hat-xi_hat-alpha_hat)^2)
  
  return(list(beta=beta_hat,active_ind=active_ind,alpha=alpha_hat,theta=theta_hat,r=r,
              sigma2=sigma2_hat,act_cov_hat=cov_hat,act.coef.mat=coef.mat,
              fitted=fitted.Y,chisq.val=chisq.val,chisq.p=chisq.p,lambda_opt=lambda_opt,
              project_measures=project_measures,out_mse=min_mse,in_mse=in_mse))
  
}






