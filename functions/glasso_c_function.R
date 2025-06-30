## Function glasso with built-in constant c
### Compute the sparse inverse of a covariance matrix
### while penalizing all the entries to $c$

library(glmnet)
source("functions/est_W.R")

offdiag <- function(S) {
    diag(S) <- 0
    return(S)
}

## Input
### Covariance matrix S
### Penalizing parameter lambda
### The constant to penalized to c.  Default c=1/d*largest eigenvalue of S
### If too large, might result in slippage at the boundary of invertibility.)
### The stopping threshold thres.  Default thres=1e-4.
### Will stop if the average absolute update of the parameter matrix (W)
### is less than thres * mean(abs(offdiag(S)))
### The maximum number of iteration.  Default iter.max=1000.

glasso_c <- function(S, lambda, c = 0, thres = 1e-2,
 iter.max = 1000) {

  tol <- thres * mean(abs(offdiag(S)))
  d <- nrow(S)
  ## Start the clock!
  # ptm <- proc.time()

  ## Initial guess
  W <- S 

  ## Iteration
  delta <- 1
  i <- 1
  while ((delta > tol) && i <= iter.max) {
    i <- i + 1
    W_old <- W
    for (j in 1:d) {
      W11 <- W[-j, -j]
      w12 <- W[j, -j]
      w22 <- W[j, j]
      w22star <- S[j, j]
      s12 <- S[j, -j]
      Theta11Inv <- W11 - w12 %*% t(w12) / w22
      eigen.decomp <- eigen(Theta11Inv)

      # diag.value <- pmax(eigen.decomp$values, 0.0001)
      diag.value <- eigen.decomp$values
      # print(diag.value)
      Q <- eigen.decomp$vectors
      A <- Q %*% diag(sqrt(diag.value)) %*% t(Q)
      AInv <- Q %*% diag(sqrt(diag.value)^(-1)) %*% t(Q)
      b <- AInv %*% (- s12 - Theta11Inv %*% rep(c * w22star, d - 1))
      lasso.obj <- glmnet(A, b, alpha = 1, lambda = lambda / (d-1),
        standardize = F, intercept = F)
      alpha <- as.vector(coef(lasso.obj))[-1]
      # rm(lasso.obj)
      # gc()
      theta12 <- alpha / w22star + rep(c, d - 1)
      # theta22 <- 1 / w22star + t(theta12) %*% Theta11Inv %*% theta12
      tmp_vec <- Theta11Inv %*% theta12
      W11 <- Theta11Inv + w22star * tmp_vec %*% t(tmp_vec)
      w12 <- - w22star * tmp_vec
      w22 <- w22star
      
      W[-j, -j] <- W11
      W[j, -j] <- w12
      W[-j, j] <- w12
      W[j, j] <- w22
      W <- (W + t(W)) / 2
    }
    
    delta <- mean(abs(W - W_old))
  }
  
  Theta_hat <- solve(W) 

  return(list(Theta = Theta_hat, Sigma = W))
}

glasso_c_reest <- function(S, lambda, c = 0, thres = 1e-2,
  iter.max = 1000) {

  if (c == 0) c <- 1 / (d * eigen(S)$values[d-1])

  d <- nrow(S)


  ## Initial guess
  S <- S + 1 / c / (d ^ 2) * matrix (1, d, d)

  # Theta_simple <- solve(S) - c * matrix (1, d, d)
  # Theta_demean <- demean_matrix(Theta_simple)
  # Theta_demean <- Theta_simple

  res <- glasso_c(S, lambda, c, 1e-2, iter.max)

  # W <- res$Sigma
  # W <- W - 1 / c / (d ^ 2) * matrix (1, d, d)

  # Theta_hat <- trim_matrix(res$Theta, c)
  Theta_hat <- res$Theta - c * matrix (1, d, d)

  select <- abs(res$Theta - c * matrix (1, d, d)) < thres
  # select <- abs(Theta_hat) < thres
  for (i in 1:d) select[i, i] <- FALSE

  # Theta_demean[select] <- 0
  # Theta_hat[select] <- 0

  return(list(Theta_hat = Theta_hat, graph = select))

}