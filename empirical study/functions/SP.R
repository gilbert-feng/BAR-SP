# intermediate  function to implement SP estimation (free of repeated SVD)
## mainly based on the functions provided in R package "NetworkReg"
SP.Prod <- function(X,x_proj,Z_hat,Y,W_cup,K,r,homo,alpha.CI) {
  # X: covariate matrix
  # x_proj: covariate-related projection matrix
  # Z_hat: rotated orthonormal covariate matrix
  # Y: response variable
  # W_cup: rotated network spectral matrix
  # K: the dimension of network spectral information
  # r: dimension of common component
  # homo: 1: homoskedastic disturbance 0: heteroskedastic disturbance
  # alpha.CI: significance level


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
    xi_hat <- 0
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

  
  
  # BAR estimate with optimal lambda
  beta_hat <- x_proj%*%beta_z
  active_ind <- 1:ncol(X)
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
  est_error <- Y-X%*%beta_hat-xi_hat-alpha_hat
  
  in_mse <- mean(est_error^2)


  return(list(beta=beta_hat,active_ind=active_ind,alpha=alpha_hat,theta=theta_hat,r=r,
              sigma2=sigma2_hat,act_cov_hat=cov_hat,act.coef.mat=coef.mat,
              fitted=fitted.Y,chisq.val=chisq.val,chisq.p=chisq.p,project_measures=project_measures,est_error=est_error))

}

# intermediate function to implement SP estimation, but provide additional prediction measures (mainly used in model prediction comparison simulations)
SP.Prod.pred <- function(X,x_proj,Z_hat,Y,W_cup,K,r,homo,K_fold=2,alpha.CI, regular = 0) {
  # regular: whether add small perturbation to make OLS stable
  
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
    xi_hat <- 0
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
  
  pred_mse <- rep(0,K_fold)
  #out_mse <- rep(0,K_fold)
  groups <- sample(rep(1:K_fold, length.out = n))
  
  for (fold_id in 1:K_fold) {
    test_index <- (groups == fold_id)
    train_index <- !test_index
    X_train <- X[train_index, ]
    X_test <- X[test_index, ]
    
    # Training pseudo response
    PC_y_train <- P_C[train_index, train_index] %*% Y[train_index]
    beta_new <- solve(crossprod(X_train) + 1e-9*diag(ncol(X_train))*regular,t(X_train))%*%PC_y_train
    
    # Prediction MSE
    y_pred <- X_test %*% beta_new
    y_proj <- P_C[test_index, test_index] %*% Y[test_index]
    pred_mse[fold_id] <- mean((y_proj - y_pred)^2)
    #out_mse[fold_id] <- mean((Y[test_index]-y_pred-xi_hat[test_index]-alpha_hat[test_index])^2)
  
  }
  min_mse <- min(pred_mse)
  #out_min_mse <- min(out_mse)
  
  # BAR estimate with optimal lambda
  beta_hat <- x_proj%*%beta_z
  active_ind <- 1:ncol(X)
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
  
  est_error <- Y-X%*%beta_hat-xi_hat-alpha_hat
  
  in_mse <- mean(est_error^2)
  
  return(list(beta=beta_hat,active_ind=active_ind,alpha=alpha_hat,theta=theta_hat,r=r,
              sigma2=sigma2_hat,act_cov_hat=cov_hat,act.coef.mat=coef.mat,
              fitted=fitted.Y,chisq.val=chisq.val,chisq.p=chisq.p,
              project_measures=project_measures,out_mse=min_mse,in_mse=in_mse,est_error=est_error))
  
}


