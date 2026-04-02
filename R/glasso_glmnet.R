# Graphical LASSO with glmnet as the inner solver
# Outer loop follows the same column-cycling algorithm as glasso_c_cpp,
# but each column's lasso subproblem is solved via glmnet::glmnet.

library(glmnet)

#' Graphical LASSO using glmnet for inner lasso
#'
#' @param S Covariance matrix (d x d)
#' @param lambda Penalization parameter
#' @param c Constant to penalize towards (default 0)
#' @param iter.max Maximum number of outer iterations (default 200)
#' @return List with Theta, Sigma, graph (logical matrix, TRUE = edge absent)
glasso_c_glmnet <- function(S, lambda, c = 0, iter.max = 200) {
  d <- nrow(S)

  # Convergence tolerance from scale of S
  S_offdiag <- S
  diag(S_offdiag) <- 0
  tol <- 1e-2 * mean(abs(S_offdiag))

  W <- S
  graph <- matrix(0L, d, d)
  Theta_direct <- matrix(0, d, d)

  for (iteration in seq_len(iter.max)) {
    W_old <- W

    for (j in seq_len(d)) {
      notj <- setdiff(seq_len(d), j)

      W11 <- W[notj, notj]
      w12 <- W[notj, j]
      w22 <- W[j, j]
      w22star <- S[j, j]
      s12 <- S[notj, j]

      Theta11Inv <- W11 - tcrossprod(w12) / w22

      # Linear term: h = s12 + c * w22star * Theta11Inv %*% 1
      h <- s12 + c * w22star * rowSums(Theta11Inv)

      # Construct design matrix X and response y for glmnet
      # We need: X^T X / n = Theta11Inv and X^T y / n = -h
      # where n = d - 1
      n_obs <- d - 1
      R <- tryCatch(chol(Theta11Inv), error = function(e) NULL)
      if (is.null(R)) {
        # Theta11Inv not PD; add progressively larger ridge until it works
        for (ridge in c(1e-8, 1e-6, 1e-4, 1e-2)) {
          R <- tryCatch(chol(Theta11Inv + ridge * diag(n_obs)), error = function(e) NULL)
          if (!is.null(R)) break
        }
        if (is.null(R)) {
          # Skip this column update entirely
          graph[notj, j] <- 1L
          next
        }
      }
      # R is upper triangular: t(R) %*% R = Theta11Inv
      X <- sqrt(n_obs) * t(R)  # n_obs x n_obs, lower triangular
      y <- -sqrt(n_obs) * forwardsolve(t(R), h)

      # Solve with glmnet (single lambda, no standardization, no intercept)
      fit <- glmnet(X, y, lambda = lambda, standardize = FALSE,
                    intercept = FALSE, thresh = 1e-7)
      alpha <- as.vector(fit$beta)

      theta12 <- alpha / w22star + c

      # Update W (same formulas as C++ version)
      tmp_vec <- as.vector(Theta11Inv %*% theta12)
      W[notj, notj] <- Theta11Inv + w22star * tcrossprod(tmp_vec)
      W[notj, j] <- -w22star * tmp_vec
      W[j, notj] <- -w22star * tmp_vec
      W[j, j] <- w22star
      W <- (W + t(W)) / 2

      # Track sparsity
      graph[notj, j] <- as.integer(alpha == 0)

      # Store column solution
      Theta_direct[notj, j] <- theta12
      Theta_direct[j, j] <- 1 / w22star + sum(theta12 * tmp_vec)
    }

    delta <- mean(abs(W - W_old))
    if (delta < tol) break
  }

  Theta <- tryCatch(solve(W), error = function(e) {
    (Theta_direct + t(Theta_direct)) / 2
  })

  # Symmetrize graph
  graph <- graph * t(graph)
  diag(graph) <- 0L
  graph <- graph == 1L

  list(Theta = Theta, Sigma = W, graph = graph)
}

#' Graphical LASSO with re-estimation using glmnet
#'
#' @param S Covariance matrix (d x d)
#' @param lambda Penalization parameter
#' @param c Constant (default 0, auto-computed from eigenvalues)
#' @param iter.max Maximum iterations (default 200)
#' @return List with Theta_hat (precision matrix) and graph (logical matrix)
glasso_c_reest_glmnet <- function(S, lambda, c = 0, iter.max = 200) {
  d <- nrow(S)

  if (c == 0) {
    eigvals <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    pos_eigvals <- eigvals[eigvals > 1e-10]
    c <- 1 / (d * min(pos_eigvals))  # smallest positive eigenvalue
  }

  # Modify S
  S_mod <- S + (1 / c / d^2) * matrix(1, d, d)

  res <- glasso_c_glmnet(S_mod, lambda, c, iter.max)
  Theta_hat <- res$Theta - c
  select <- res$graph

  # Zero out entries identified as absent
  Theta_hat[select] <- 0

  list(Theta_hat = Theta_hat, graph = select)
}
