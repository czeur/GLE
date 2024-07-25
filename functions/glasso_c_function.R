## Function glasso with built-in constant c
### Compute the sparse inverse of a covariance matrix while penalizing all the entries to $c$
library(glmnet)

offdiag <- function(S){
    diag(S) <- 0
    return(S)
}
## Input
### Covariance matrix S 
### Penalizing parameter lambda
### The constant to penalized to c.  Default c=1/d*largest eigenvalue of S.  If too large, might result in slippage at the boundary of invertibility.)
### The stopping threshold thres.  Default thres=1e-4.  Will stop if the average absolute update of the parameter matrix (W) < thres * mean(abs(offdiag(S)))
### The maximum number of iteration.  Default iter.max=1000.

glasso_c <- function(S,lambda,c=1/(d*eigen(S)$values[1]),thres=1e-4,iter.max=1000){

  d <- nrow(S)

  tol <- thres * mean(abs(offdiag(S)))
  
  ## Start the clock!
  ptm <- proc.time()
  
  ## Initial guess
  W <- S
  
  ## Iteration
  delta <- 1
  i <- 1
  while ((delta>tol) & i<=iter.max ){
    i <- i+1
    W_old <- W
    for (j in 1:d){
      W11 <- W[-j,-j]
      w12 <- W[j,-j]
      w22 <- W[j,j]
      w22star <- S[j,j]
      s12 <- S[j,-j]
      Theta11Inv <- W11 - w12%*%t(w12)/w22
      eigen.decomp <- eigen(Theta11Inv)
      Q <- eigen.decomp$vectors
      A <- Q%*%diag(sqrt(eigen.decomp$values))%*%t(Q)
      AInv <- Q%*%diag(sqrt(eigen.decomp$values)^(-1))%*%t(Q)
      b <- AInv%*% (- s12 - Theta11Inv%*%rep(c*w22star,d-1))
      lasso.obj <- glmnet(A, b, alpha = 1, lambda = lambda/(d-1), standardize = F, intercept = F)
      alpha <- as.vector(coef(lasso.obj))[-1]
      theta12 <- alpha/w22star + rep(c,d-1)
      theta22 <- 1/w22star + t(theta12)%*%Theta11Inv%*%theta12
      tmp_vec <- Theta11Inv %*% theta12
      W11 <- Theta11Inv + w22star * tmp_vec %*% t(tmp_vec)
      w12 <- - w22star * tmp_vec
      w22 <- w22star
      
      W[-j,-j] <- W11
      W[j,-j] <- w12
      W[-j,j] <- w12
      W[j,j] <- w22
      W <- (W + t(W))/2
    }
    
    print(i)
    delta <- mean(abs(W - W_old))
  }
  
  ## The timer
  print(proc.time() - ptm)
  Theta <- solve(W)
  return(list(Theta=Theta,Sigma=W))
}