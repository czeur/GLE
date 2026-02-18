# Incoherence measures for graphical models

library(graphicalExtremes)
library(igraph)

incoherence_pess <- function(Sig, S, Sc) {
  if (length(Sc) != 0)
    1 - max(sapply(Sc, function(e) sum(abs(Sig[e, S] %*% solve(Sig[S, S])))))
  else
    1
}

incoherence <- function(Sig, S, Sc, signs) {
  if (length(Sc) != 0)
    1 - max(abs(sapply(Sc, function(e)
      sum(Sig[e, S] %*% solve(Sig[S, S]) %*% signs))))
  else
    1
}

#' Compute GL and NS incoherence for a precision matrix
#'
#' @param Sig Covariance (or correlation) matrix
#' @param cor.scale Scale to correlation (default TRUE)
#' @param pessimistic Use pessimistic NS bound (default FALSE)
#' @return Length-2 vector: c(GLi, NSi)
GLNSi <- function(Sig, cor.scale = TRUE, pessimistic = FALSE) {
  d <- ncol(Sig)
  if (cor.scale) {
    D_inv <- diag(Sig)^(-0.5) * diag(d)
    Sig <- D_inv %*% Sig %*% D_inv
  }

  Th <- solve(Sig)
  Th <- Th * (abs(Th) >= 1e-6)

  gr <- graph_from_adjacency_matrix(Th != 0, mode = "undirected")
  S <- which(as.vector(gr[1:d]) == 1)
  Sc <- which(as.vector(gr[1:d]) == 0)
  H <- kronecker(Sig, Sig)
  GLi <- incoherence_pess(H, S, Sc)

  NSi <- min(sapply(1:d, function(a) {
    if (sum(Th[-a, a] != 0) == 0) return(1)
    S <- which(Th[, a] != 0)
    S <- S[S != a]
    Sc <- which(Th[, a] == 0)
    if (pessimistic) incoherence_pess(Sig, S, Sc)
    else incoherence(Sig, S, Sc, sign(Th[S, a]))
  }))

  c(GLi, NSi)
}

#' Compute incoherence from a variogram matrix
#'
#' @param G Variogram matrix
#' @param cor.scale Scale to correlation (default TRUE)
#' @param pessimistic Use pessimistic bound (default FALSE)
#' @return List with GLi and NSi vectors (length d)
Gamma2Inc <- function(G, cor.scale = TRUE, pessimistic = FALSE) {
  Inc <- t(sapply(1:nrow(G), function(m) {
    Sig <- Gamma2Sigma(G, m)
    GLNSi(Sig, cor.scale, pessimistic)
  }))
  list(GLi = Inc[, 1], NSi = Inc[, 2])
}
