# R wrapper for C++ W matrix estimation

Rcpp::sourceCpp("src/est_W.cpp")

#' Estimate the W matrix from data via rank transformation + threshold truncation
#'
#' @param data Data matrix (n x d)
#' @param q.threshold Quantile threshold (e.g. 0.9)
#' @return List with cov (trimmed covariance matrix) and subcovlist
#'   (list of d conditional covariance matrices)
est_W <- function(data, q.threshold) {
  est_W_cpp(data, q.threshold)
}
