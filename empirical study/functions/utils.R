
# function to calculate the precision and recall
ZeroCount <- function(beta_0, beta_hat, threshold) {
  beta_hat_index <- ifelse(abs(beta_hat) < threshold, 0, 1)
  beta_0_index <- ifelse(abs(beta_0) < threshold, 0, 1)
  
  # Confusion elements
  TP <- sum(beta_0_index == 0 & beta_hat_index == 0)  
  FP <- sum(beta_0_index == 1 & beta_hat_index == 0)  
  FN <- sum(beta_0_index == 0 & beta_hat_index == 1)  
  
  # Smoothing
  Precision <- (TP + 1) / (TP + FP + 1)
  Recall <- (TP + 1) / (TP + FN + 1)
  
  return(list(Precision = Precision, Recall = Recall))
}

# function to calculate the G index
GroupCount <- function(beta_hat, G, r, threshold) {
  p <- length(beta_hat)
  beta_hat_index <- rep(1, p)
  beta_hat_index[abs(beta_hat[1:p]) < threshold] <- 0
  G_num <- (p-r)/G
  G_total <- numeric(G)
  
  for (i in 1:G) {
    temp_ind <- ((i-1)*G_num+1+r):((i-1)*G_num+1+r+G_num-1)
    if (i<=G/2) {
      G_total[i] <- as.integer(all(beta_hat_index[temp_ind] != 0))
    } else {
      G_total[i] <- as.integer(all(beta_hat_index[temp_ind] == 0))
      
    }
  }

  G_measure <- mean(G_total)
  
  return(G_measure)
}

# function to generate network based on the linking probability matrix P
## collected from the R package "NetworkReg"
net.gen.from.P <- function(P,mode="undirected"){
  n <- nrow(P)
  if(mode=="undirected"){
    upper.index <- which(upper.tri(P))
    upper.p <- P[upper.index]
    upper.u <- runif(n=length(upper.p))
    upper.A <- rep(0,length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    
    
    A <- matrix(0,n,n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    
  }else{
    A <- matrix(0,n,n)
    R <- matrix(runif(n^2),n,n)
    A[R<P] <- 1
    diag(A) <- 0
  }
  return(A)
}

# function to estimate the SBM network
## collected from the R package "randnet"
SBM.estimate <- function(A,g){
  n <- nrow(A)
  K <- length(unique(g))
  B <- matrix(0,K,K)
  M <- matrix(0,n,K)
  for(i in 1:K){
    for(j in i:K){
      if(i!=j){
        B[j,i] <- B[i,j] <- mean(A[which(g==i),which(g==j)])
      }else{
        n.i <- length(which(g==i))
        if(n.i>1){
          B[i,i] <- sum(A[which(g==i),which(g==i)])/(n.i^2 - n.i)
        }else{
          B[i,i] <- 0
        }
      }
    }
  }
  M[matrix(c(1:n,g),ncol=2)] <- 1
  P <- M%*%B%*%t(M)
  return(list(B=B,Phat=P,g=g))
}

# function to generate SBM network
## collected from the R package "randnet"
BlockModel.Gen <- function(lambda,n,beta=0,K=3,w=rep(1,K),Pi=rep(1,K)/K,rho=0,simple=TRUE,power=TRUE,alpha=5,degree.seed=NULL){
  P0 <- diag(w)
  if(beta > 0){
    P0 <- matrix(1,K,K)
    diag(P0) <- w/beta
  }
  Pi.vec <- matrix(Pi,ncol=1)
  P <- P0
  
  M <- matrix(0,n,K)
  membership <- sample(x=K,size=n,replace=TRUE,prob=Pi)
  M[cbind(1:n,membership)] <- 1
  A.bar <- M%*%P%*%t(M)
  node.degree <- rep(1,n)
  if(rho>0){
    if(simple){
      node.degree[runif(n)<rho] <- 0.2
    }else{
      if(power==FALSE){
        node.degree <- runif(n)*0.8 + 0.2
      }else{
        MM <- ceiling(n/300)
        if(is.null(degree.seed)){
          degree.seed <- rplcon(300,1,alpha)
        }### sample fixed number of values from power law and then randomly assign to be degrees. Easier to control noises across different sample sizes.
        node.degree <- sample(degree.seed,size=n,replace=TRUE)
      }
    }}
  A.bar <- t(t(A.bar*node.degree)*node.degree)
  A.bar <- A.bar*lambda/mean(colSums(A.bar))
  #diag(A.bar) <- 0
  #avg.d <- mean(colSums(A.bar))
  #A.bar <- A.bar*lambda/avg.d
  upper.index <- which(upper.tri(A.bar))
  upper.p <- A.bar[upper.index]
  upper.u <- runif(n=length(upper.p))
  upper.A <- rep(0,length(upper.p))
  upper.A[upper.u < upper.p] <- 1
  A <- matrix(0,n,n)
  A[upper.index] <- upper.A
  A <- A + t(A)
  diag(A) <- 0
  return(list(A=A,g=membership,P=A.bar,theta=node.degree))
}

# function to determine the dimension of common component and calculate
# the corresponding projection operator
ProdSVD <- function(X,A,K,r=NULL,thr=NULL,boot.thr=TRUE,boot.n=50) {
  n <- nrow(X)
  p <- ncol(X)
  eigen_A <- eigs_sym(A=A,k=K)
  X_svd <- svd(X)
  x_proj <- X_svd$v%*%(t(X_svd$u)/X_svd$d)
  Z <- X_svd$u
  W <- matrix(eigen_A$vectors[,1:K],ncol=K)
  ProdSVD <- svd(t(Z)%*%W,nu=p,nv=K)
  S_hat <- ProdSVD$d
  Z_hat <- Z%*%ProdSVD$u
  W_cup <- W%*%ProdSVD$v
  
  # determine the dimension of common component
  if(is.null(r)){
    if(is.null(thr)){
      if(boot.thr){
        P0 <- eigen_A$vectors %*%t(eigen_A$vectors*eigen_A$values)
        P0[P0>1] <- 1
        P0[P0<0] <- 0
        P0 <- P0/mean(rowSums(P0))*mean(rowSums(A))
        eig_P0 <- eigs_sym(P0,k=K)
        fake.svd <- svd(t(Z)%*%matrix(eig_P0$vectors[,1:K],ncol=K))
        message("Use bootstrapping to find r-threshold....")
        sim.s <- matrix(0,boot.n,min(p,K))
        for(II in 1:boot.n){
          A.sim <- net.gen.from.P(P0)
          eig.A.sim <- eigs_sym(A.sim,k=K)
          sim.svd <- svd(t(Z)%*%matrix(eig.A.sim$vectors[,1:K],ncol=K))
          sim.s[II,] <- sim.svd$d
        }
        
        thr <- 1-max(abs(t(sim.s)-fake.svd$d))
        message(paste("Select r by threshold",thr))
        
      }else{
        dhat <- max(rowSums(A))
        thr <- 1 - 4*sqrt(p*K*log(n))/dhat
        if(thr < 0.05){
          thr <- 0.05
        }
        message(paste("Select r by asymptotic threshold",thr))
      }
    }
    r <- sum(ProdSVD$d>=thr)
    message(paste("Select r =",r))
  }
  
  return(list(x_proj = x_proj, S_hat = S_hat,
              Z = Z, Z_hat = Z_hat,
              W = W, W_cup = W_cup,
              r = r, thr = thr))
}

# function to implement the broken adaptive ridge estimation
BAR.iter <- function(X,y,lambda,K_fold=2,iter_max=1e3,conv_crit=1e-5, regular = 0) {
  
  XX <- crossprod(X)
  n <- nrow(X)
  p <- ncol(X)
  # Initial estimator
  beta_old <- solve(XX + 1e-9*diag(p)*regular, t(X) %*% y)

  iter_num <- 1
  converged <- FALSE
  
  # Iteration
  while (iter_num < iter_max && !converged) {
    # Calculate truncated diagonal weight
    diag_vec <- pmin(1e9, lambda / (beta_old^2))
    XX_new <- XX + diag(diag_vec, nrow = ncol(X))
    # New estimator
    beta_new <- solve(XX_new, t(X) %*% y)
    
    # Check convergence
    if (sum(abs(beta_new - beta_old)) < conv_crit) {
      converged <- TRUE
    } else {
      beta_old <- beta_new
      iter_num <- iter_num + 1
    }
  }
  return(beta_new)
}

