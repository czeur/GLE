# R wrappers for C++ graphical LASSO implementation

Rcpp::sourceCpp("src/glasso_c.cpp")

#' Graphical LASSO with built-in constant c
#'
#' Compute the sparse inverse of a covariance matrix while penalizing
#' all off-diagonal entries towards c.
#'
#' @param S Covariance matrix (d x d)
#' @param lambda Penalization parameter
#' @param c Constant to penalize towards (default 0)
#' @param iter.max Maximum number of iterations (default 1000)
#' @return List with Theta (precision matrix), Sigma (covariance estimate),
#'   and graph (logical matrix, TRUE = edge absent)
glasso_c <- function(S, lambda, c = 0, iter.max = 1000) {
  res <- glasso_c_cpp(S, lambda, c, iter.max)
  res$graph <- res$graph == 1
  res
}

#' Graphical LASSO with re-estimation and sparsity detection
#'
#' Runs glasso_c with the constant c modification, then identifies the
#' sparsity pattern (zero entries) from the LASSO exact zeros.
#'
#' @param S Covariance matrix (d x d)
#' @param lambda Penalization parameter
#' @param c Constant (default 0, auto-computed from eigenvalues)
#' @param iter.max Maximum iterations (default 1000)
#' @return List with Theta_hat (precision matrix) and graph (logical matrix,
#'   TRUE = edge absent)
glasso_c_reest <- function(S, lambda, c = 0, iter.max = 1000) {
  res <- glasso_c_reest_cpp(S, lambda, c, iter.max)
  # C++ returns graph as numeric 0/1; convert to logical for safe R indexing
  res$graph <- res$graph == 1
  res
}
