# Extended eglearn and neighborhood selection functions
# Ported from applications/functions_paper.R

library(Matrix)
library(graphicalExtremes)
library(igraph)

mychol <- function(M) {
  d <- nrow(M)
  n <- rankMatrix(M)
  if (n == d) return(chol(M))
  R <- chol(M, pivot = TRUE)
  o <- order(attr(R, "pivot"))
  R[1:n, o]
}

#' Neighborhood selection via LASSO (Meinshausen-Buhlmann)
#'
#' @param data Data matrix (from Cholesky decomposition of Sigma_k)
#' @param samp_size Effective sample size
#' @param lambda Vector of penalty parameters
#' @param thr_zero Threshold for zeroing coefficients
#' @param ic Whether to compute information criteria
#' @param refit Whether to refit without penalty for IC computation
#' @return List with adj.est (and adj.ic.est if ic=TRUE)
glasso_mb2 <- function(data, samp_size, lambda, thr_zero = 1e-10,
                       ic = FALSE, refit = TRUE) {
  dd <- ncol(data)
  S_tmp <- t(data) %*% data
  data_std <- data %*% diag(diag(S_tmp)^(-1/2))
  adj.est <- array(NA, dim = c(dd, dd, length(lambda)))
  if (ic) adj.ic.est <- array(NA, dim = c(dd, dd, 3))
  lambda_order <- order(lambda, decreasing = TRUE)
  lambda_dec <- sort(lambda, decreasing = TRUE)

  for (i in 1:dd) {
    X <- data_std[, -i]
    Y <- data_std[, i]
    lasso_fit <- glmnet::glmnet(
      x = X, y = Y, family = "gaussian",
      lambda = lambda_dec / nrow(X) * samp_size / (samp_size - 1),
      standardize = FALSE, intercept = FALSE)

    if (i == 1) {
      null.vote <- array(0, dim = c(dd, dd, length(lambda)))
      if (ic) null.vote.ic <- array(0, dim = c(dd, dd, 3))
    }

    null.vote[i, -i, ] <- null.vote[i, -i, ] +
      (abs(as.matrix(lasso_fit$beta)) <= thr_zero)
    null.vote[-i, i, ] <- null.vote[-i, i, ] +
      (abs(as.matrix(lasso_fit$beta)) <= thr_zero)

    if (ic) {
      if (refit) {
        dev <- sapply(1:length(lambda), function(l) {
          to_excl <- which(abs(as.matrix(lasso_fit$beta)[, l]) <= 1e-10)
          if (length(to_excl) == ncol(X)) return(samp_size * sum(Y^2))
          X_tmp <- as.matrix(X[, -to_excl])
          samp_size * sum(((diag(nrow(X_tmp)) -
            X_tmp %*% solve(t(X_tmp) %*% X_tmp) %*% t(X_tmp)) %*% Y)^2)
        })
      } else {
        dev <- samp_size * (1 - lasso_fit$dev.ratio) * lasso_fit$nulldev
      }

      dfs <- lasso_fit$df
      aic_fn <- function(n, p) 2
      bic_fn <- function(n, p) log(n)
      mbic_fn <- function(n, p) log(n) * log(log(p))

      aic.idx <- which.min(dev + aic_fn(samp_size, dd) * dfs)
      bic.idx <- which.min(dev + bic_fn(samp_size, dd) * dfs)
      mbic.idx <- which.min(dev + mbic_fn(samp_size, dd) * dfs)

      null.vote.ic[i, -i, ] <- null.vote.ic[i, -i, ] +
        (abs(as.matrix(lasso_fit$beta[, c(aic.idx, bic.idx, mbic.idx)])) <= thr_zero)
      null.vote.ic[-i, i, ] <- null.vote.ic[-i, i, ] +
        (abs(as.matrix(lasso_fit$beta[, c(aic.idx, bic.idx, mbic.idx)])) <= thr_zero)
    }
  }

  adj.est[, , lambda_order] <- null.vote <= 1
  if (ic) {
    adj.ic.est <- null.vote.ic <= 1
    list(adj.est = adj.est, adj.ic.est = adj.ic.est)
  } else {
    list(adj.est = adj.est)
  }
}

#' Extended eglearn with HR likelihood-based model selection
#'
#' @param data Data matrix
#' @param p Threshold probability (NULL if data already mpareto)
#' @param rholist Vector of regularization parameters
#' @param thr_zero Threshold for zero detection
#' @param reg_method "ns" (neighborhood selection) or "glasso"
#' @param ic "hr" (HR likelihood) or "ns" (Gaussian IC per regression)
#' @param refit Refit without penalty for IC
#' @param return_Gamma Whether to return estimated Gamma matrices
#' @return List with graph, Gamma, rholist, graph_ic, Gamma_ic
eglearn2 <- function(data, p = NULL,
                     rholist = c(0.1, 0.15, 0.19, 0.205),
                     thr_zero = 1e-5,
                     reg_method = c("ns", "glasso"),
                     ic = c("hr", "ns"),
                     refit = FALSE,
                     return_Gamma = FALSE) {
  reg_method <- match.arg(reg_method)
  ic <- match.arg(ic)
  if (ic == "hr") return_Gamma <- TRUE

  if (any(rholist < 0)) stop("rholist must be non-negative.")

  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  Gamma <- emp_vario(data = data.std)
  sel_methods <- c("aic", "bic", "mbic")

  r <- length(rholist)
  n <- nrow(data.std)
  d <- ncol(Gamma)

  null.vote <- array(0, dim = c(d, d, r))
  null.vote.ic <- array(0, dim = c(d, d, length(sel_methods)))

  for (k in 1:d) {
    Sk <- Gamma2Sigma(Gamma = Gamma, k = k)
    if (reg_method == "glasso") {
      gl.fit <- lapply(seq_along(rholist), function(i) {
        glassoFast::glassoFast(S = Sk, rho = rholist[i], thr = 1e-8,
                               maxIt = 100000)$wi
      })
      gl.tmp <- array(unlist(gl.fit), dim = c(d - 1, d - 1, r))
      null.vote[-k, -k, ] <- null.vote[-k, -k, , drop = FALSE] +
        (abs(gl.tmp) <= thr_zero)
    } else if (reg_method == "ns") {
      samp_size <- length(which(data.std[, k] > 1))
      X <- mychol(Sk)
      gl.tmp <- glasso_mb2(data = X, samp_size = samp_size, lambda = rholist,
                           thr_zero = thr_zero, ic = (ic == "ns"), refit = refit)
      null.vote[-k, -k, ] <- null.vote[-k, -k, , drop = FALSE] + (!gl.tmp$adj.est)
      if (ic == "ns") {
        null.vote.ic[-k, -k, ] <- null.vote.ic[-k, -k, , drop = FALSE] +
          (!gl.tmp$adj.ic.est)
      }
    }
  }

  adj.est <- (null.vote / (d - 2)) < 0.5
  adj.ic.est <- (null.vote.ic / (d - 2)) < 0.5

  graphs <- list()
  Gammas <- list()
  rhos <- list()
  if (ic == "hr") {
    logliks <- rep(-Inf, r)
    n_edges <- rep(Inf, r)
  }

  for (j in 1:r) {
    rho <- rholist[j]
    est_graph <- graph_from_adjacency_matrix(
      adj.est[, , j], mode = "undirected", diag = FALSE)

    Gamma_curr <- NA
    if (return_Gamma & is_connected(est_graph)) {
      try(Gamma_curr <- complete_Gamma(graph = est_graph, Gamma = Gamma),
          silent = TRUE)
    }

    if (ic == "hr" & is_connected(est_graph)) {
      if (is_valid_Gamma(Gamma_curr)) {
        logliks[j] <- graphicalExtremes:::logLH_HR(data = data.std, Gamma = Gamma_curr)
        n_edges[j] <- length(E(est_graph))
      }
    }

    graphs[[j]] <- est_graph
    Gammas[[j]] <- Gamma_curr
    rhos[[j]] <- rho
  }

  graphs_ic <- list(aic = NA, bic = NA, mbic = NA)
  Gammas_ic <- list(aic = NA, bic = NA, mbic = NA)

  if (reg_method == "ns" & ic == "ns") {
    for (l in seq_along(sel_methods)) {
      est_graph <- graph_from_adjacency_matrix(
        adj.ic.est[, , l], mode = "undirected", diag = FALSE)
      Gamma_curr <- NA
      if (return_Gamma) {
        Gamma_curr <- complete_Gamma(graph = est_graph, Gamma = Gamma)
      }
      graphs_ic[[l]] <- est_graph
      Gammas_ic[[l]] <- Gamma_curr
    }
    graph_ic <- list(aic = graphs_ic[[1]], bic = graphs_ic[[2]], mbic = graphs_ic[[3]])
    Gamma_ic <- list(aic = Gammas_ic[[1]], bic = Gammas_ic[[2]], mbic = Gammas_ic[[3]])
  } else if (ic == "hr" & any(is.finite(logliks))) {
    aic_id <- which.max(logliks - 2 * n_edges)
    bic_id <- which.max(logliks - log(n) * n_edges)
    mbic_id <- which.max(logliks - log(n) * max(1, log(log(d))) * n_edges)
    graph_ic <- list(aic = graphs[[aic_id]], bic = graphs[[bic_id]], mbic = graphs[[mbic_id]])
    Gamma_ic <- list(aic = Gammas[[aic_id]], bic = Gammas[[bic_id]], mbic = Gammas[[mbic_id]])
  } else {
    graph_ic <- list(aic = NULL, bic = NULL, mbic = NULL)
    Gamma_ic <- list(aic = NULL, bic = NULL, mbic = NULL)
  }

  list(graph = graphs, Gamma = Gammas, rholist = rhos,
       graph_ic = graph_ic, Gamma_ic = Gamma_ic)
}
