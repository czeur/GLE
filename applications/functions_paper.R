mychol <- function(M){
  d <- nrow(M)
  n <- rankMatrix(M)
  if(n==d) return(chol(M))
  else{
    R <- chol(M, pivot=TRUE)
    o <- order(attr(R, "pivot"))
    return(R[1:n,o])
  }
}

aic <- function(n, p) 2
bic <- function(n, p) log(n)
mbic <- function(n, p) log(n) * log(log(p)) # modified BIC of Wang & Leng, JRSSB 2009



glasso_mb2 <- function(data, samp_size, lambda, thr_zero = 1e-10, ic = FALSE, refit = TRUE){
  
  # Initialize variables
  dd <- ncol(data)
  S_tmp <- t(data) %*% data
  data_std <- data %*% diag(diag(S_tmp)^(-1/2))
  adj.est <- array(NA, dim = c(dd, dd, length(lambda)))
  if(ic) adj.ic.est <- array(NA, dim = c(dd, dd, 3))
  lambda_order <- order(lambda, decreasing = TRUE)
  lambda_dec <- sort(lambda, decreasing = TRUE)
  
  # Loop through variables
  for(i in (1:dd)){
    X <- data_std[,-i]
    Y <- data_std[,i]
    lasso_fit <- glmnet::glmnet(x = X, y = Y, family = "gaussian",
                                lambda = lambda_dec/nrow(X)*samp_size/(samp_size-1),
                                standardize = F, intercept = F)
    if(i==1){
      null.vote <- array(0, dim = c(dd, dd, length(lambda)))
      if(ic) null.vote.ic <- array(0, dim = c(dd, dd, 3))
    }
    null.vote[i, -i, ] <- null.vote[i, -i, ] +
      (abs(as.matrix(lasso_fit$beta)) <= thr_zero)
    null.vote[-i, i, ] <- null.vote[-i, i, ] +
      (abs(as.matrix(lasso_fit$beta)) <= thr_zero)
    
    if(ic){
      # refitting the estimated models without penalty to extract the likelihood value
      if(refit){
        dev <- sapply(1:length(lambda), function(l){
          to_excl <- which(abs(as.matrix(lasso_fit$beta)[,l]) <= 1e-10)
          if(length(to_excl) == ncol(X)) return(samp_size * sum(Y^2))
          else{
            # Calculating the regression residuals directly (squared norm of (I - X (X^T X)^{-1} X^T) Y)
            X_tmp <- as.matrix(X[,-to_excl])
            return(
              samp_size * sum(((diag(nrow(X_tmp)) - X_tmp %*% solve(t(X_tmp) %*% X_tmp) %*% t(X_tmp)) %*% Y)^2)
            )
          }
        })
      }
      # otherwise calculating the likelihood from the original (penalized) fit
      else dev <- samp_size * (1 - lasso_fit$dev.ratio) * lasso_fit$nulldev
      
      dfs <- lasso_fit$df
      aic.idx <- which.min( dev + aic(samp_size, dd) * dfs )
      bic.idx <- which.min( dev + bic(samp_size, dd) * dfs )
      mbic.idx <- which.min( dev + mbic(samp_size, dd) * dfs )
      
      null.vote.ic[i, -i, ] <- null.vote.ic[i, -i,] +
        (abs(as.matrix(lasso_fit$beta[,c(aic.idx, bic.idx, mbic.idx)])) <= thr_zero)
      null.vote.ic[-i, i, ] <- null.vote.ic[-i, i,] +
        (abs(as.matrix(lasso_fit$beta[,c(aic.idx, bic.idx, mbic.idx)])) <= thr_zero)
    }
  }
  
  adj.est[,,lambda_order] <- null.vote <= 1
  if(ic){
    adj.ic.est <- null.vote.ic <= 1
    return(list(adj.est = adj.est, adj.ic.est = adj.ic.est))
  }
  else return(list(adj.est = adj.est))
}



# information criteria as before: ic = "ns"
# information criteria by refitting the regression models inside neighborhood selection: ic = "ns", refit = TRUE
# information criteria by HR likelihood (always selects a connected graph): ic = "hr" (default)
eglearn2 <- function(data,
                     p = NULL,
                     rholist = c(0.1, 0.15, 0.19, 0.205),
                     thr_zero = 1e-5,
                     reg_method = c("ns", "glasso"),
                     ic = c("hr", "ns"),
                     refit = FALSE,
                     return_Gamma = FALSE) {
  
  # Check arguments
  reg_method <- match.arg(reg_method)
  ic <- match.arg(ic)
  if(ic == "hr") return_Gamma <- TRUE
  
  if (any(rholist < 0)) {
    stop("The regularization parameters in `rholist` must be non-negative.",
         call. = FALSE
    )
  }
  
  # Normalize data
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  }
  else {
    data.std <- data
  }
  
  # Initialize variables
  Gamma <- emp_vario(data=data.std)
  sel_methods <- c("aic", "bic", "mbic")
  
  r <- length(rholist)
  n <- nrow(data.std)
  d <- ncol(Gamma)
  
  # votes for absence of edges
  null.vote <- array(0, dim = c(d, d, r))
  null.vote.ic <- array(0, dim = c(d, d, length(sel_methods)))
  
  # Loop through variables
  for (k in 1:d) {
    Sk <- Gamma2Sigma(Gamma = Gamma, k=k)
    if (reg_method == "glasso") {
      gl.fit <- lapply(1:length(rholist), FUN = function(i) {
        glassoFast::glassoFast(S = Sk, rho = rholist[i], thr = 1e-8, maxIt = 100000)$wi
      })
      gl.tmp <- array(unlist(gl.fit), dim = c(d - 1, d - 1, r))
      null.vote[-k, -k, ] <- null.vote[-k, -k, , drop = FALSE] +
        (abs(gl.tmp) <= thr_zero)
    }
    else if (reg_method == "ns") {
      samp_size <- length(which(data.std[, k] > 1))
      X <- mychol(Sk)
      gl.tmp <- glasso_mb2(data = X, samp_size = samp_size, lambda = rholist, thr_zero = thr_zero,
                           ic = (ic == "ns"), refit = refit)
      null.vote[-k, -k, ] <- null.vote[-k, -k, , drop = FALSE] + (!gl.tmp$adj.est)
      if(ic == "ns") null.vote.ic[-k, -k, ] <- null.vote.ic[-k, -k, , drop = FALSE] + (!gl.tmp$adj.ic.est)
    }
  }
  
  # majority voting procedure
  adj.est <- (null.vote / (d-2)) < 0.5
  # only makes sense for "ns" (if using "glasso", returns NA's)
  adj.ic.est <- (null.vote.ic / (d-2)) < 0.5
  
  
  
  # create the list of estimated graphs and Gammas (if return_Gamma = TRUE)
  graphs <- list()
  Gammas <- list()
  rhos <- list()
  if(ic == "hr"){
    logliks <- rep(-Inf, r)
    n_edges <- rep(Inf, r)
  }
  
  for (j in 1:r) {
    rho <- rholist[j]
    est_graph <- igraph::graph_from_adjacency_matrix(
      adj.est[, ,j], mode = "undirected", diag = FALSE)
    
    Gamma_curr <- NA
    if(return_Gamma & igraph::is_connected(est_graph)){
      try(
        Gamma_curr <- graphicalExtremes::complete_Gamma(graph=est_graph, Gamma=Gamma),
        silent = TRUE
      )
    }
    
    if(ic == "hr" & igraph::is_connected(est_graph)){
      if(is_valid_Gamma(Gamma_curr)){
        logliks[j] <- graphicalExtremes:::logLH_HR(data = data.std, Gamma = Gamma_curr) # for true HR ll
        n_edges[j] <- length(E(est_graph))
      }
    }
    
    graphs[[j]] <- est_graph
    Gammas[[j]] <- Gamma_curr
    rhos[[j]] <- rho
  }
  
  # create the estimated graphs and Gammas corresponding to the information criteria (if return_Gamma = TRUE)
  graphs_ic <-  list(aic = NA, bic = NA, mbic = NA)
  Gammas_ic <-  list(aic = NA, bic = NA, mbic = NA)
  
  if (reg_method == "ns" & ic == "ns") {
    for (l in seq_along(sel_methods)){
      
      est_graph <-  igraph::graph_from_adjacency_matrix(
        adj.ic.est[, ,l], mode = "undirected", diag = FALSE)
      
      if(return_Gamma){
        Gamma_curr <- graphicalExtremes::complete_Gamma(graph=est_graph, Gamma=Gamma)
      }
      else Gamma_curr <- NA
      
      graphs_ic[[l]] <- est_graph
      Gammas_ic[[l]] <- Gamma_curr
    }
    graph_ic <- list(aic = graphs_ic[[1]], bic = graphs_ic[[2]], mbic = graphs_ic[[3]])
    Gamma_ic <- list(aic = Gammas_ic[[1]], bic = Gammas_ic[[2]], mbic = Gammas_ic[[3]])
  }
  
  else if(ic == "hr" & any(is.finite(logliks))){
    aic_id <- which.max(logliks - 2 * n_edges)
    bic_id <- which.max(logliks - log(n) * n_edges)
    mbic_id <- which.max(logliks - log(n) * max(1, log(log(d))) * n_edges)
    
    graph_ic <- list(aic = graphs[[aic_id]], bic = graphs[[bic_id]], mbic = graphs[[mbic_id]])
    Gamma_ic <- list(aic = Gammas[[aic_id]], bic = Gammas[[bic_id]], mbic = Gammas[[mbic_id]])
  }
  
  else{
    graph_ic <- list(aic = NULL, bic = NULL, mbic = NULL)
    Gamma_ic <- list(aic = NULL, bic = NULL, mbic = NULL)
  }
  
  return(list(graph = graphs,
              Gamma = Gammas,
              rholist = rhos,
              graph_ic = graph_ic,
              Gamma_ic = Gamma_ic))
}



set_graph_parameters <- function(graph) {
  # set parameters
  igraph::V(graph)$color <- grDevices::adjustcolor(col = "#4477AA", alpha.f = 0.4)
  igraph::V(graph)$frame.color <- grDevices::adjustcolor(col = "#4477AA", alpha.f = 1)
  igraph::V(graph)$label.color <- "black"
  igraph::V(graph)$size <- 15
  igraph::E(graph)$width <- 2
  igraph::E(graph)$color <- "darkgrey"

  # return graph
  return(graph)
}



sim_study <- function(d = 5, 
                          n = 100,
                          p = NULL,
                          method = c("maxstable", "mpareto"), 
                          m = 2, 
                          gen_model = c("BA", "block"),
                          alphad = 1,
                          reg_method = c("ns", "glasso", "emst", "MTP2"),
                          ic = c("hr", "ns"),
                          rhostring = "seq(0.01,5,length=100)",
                          incoh = FALSE,
                          rng = NULL){
  ## perform a simulation study to measure performance of the EMTP2 block descent algorithm.
  ##
  ## Args:
  ##    - d: dimension.
  ##    - n: number of samples.
  ##    - p: threshold probability.
  ##    - method: data generation method.
  ##    - alphad: parameter for the block model.
  ##    - m: the number of edges added in each setp of the Barabasi-Albert model (m=1 is a tree)
  ##    - reg_method: regression method used to estimate the extremal graphical structure.
  ##    - ic: information criteria used for model selection, either HR likelihood ("hr", default)
  ##      or Gaussian likelihood inside each lasso regression ("ns")
  ##    - rhostring: the list of penality parameters; must be given a string.
  ##    - incoh: If TRUE, returns incoherence values.
  ## Returns:
  ##  a tibble with F scores
  
  # check arguments
  
  rholist <-  eval(parse(text = rhostring))
  gen_model <- match.arg(gen_model)
  
  F1 <- numeric(length(rholist))
  incoh_pos <- numeric(2)
  F1_ic <- NA
  
  # set seed, if applicable
  if (!is.null(rng)){
    rng_sims <- rng[[1]]
    rngtools::setRNG(rng_sims)
  }
  
  
  if(gen_model=="BA"){
    BA_model <- generate_BA_model(d = d, m = m)
    G <- BA_model$G
    g <- BA_model$graph
  }
  if(gen_model=="block"){
    block_model <- generate_block_model(ncliques = 6, clique_size = 4, alphad=alphad) 
    g <- block_model$graph
    G <- .05*block_model$G
    d <- nrow(G)
  }
  
  T <- Gamma2Theta(G)
  nb.edges <- length(E(g))
  partial_pos <- (sum(sum(T > 1e-6)) - d)/2 / nb.edges
  partial_neg <- sum(sum(T < -1e-6))/2 / nb.edges
  
  
  if(incoh==TRUE){
    incoh_pos <- sapply(Gamma2Inc(G = G, pessimistic = TRUE), FUN = function(x) mean(x))  
  }
  else incoh_pos <- c(NA,NA)
  
  # perform simulation
  if(method=="maxstable")  X <- rmstable(n=n, d=d, model="HR", par=G)
  if(method=="mpareto")  X <- rmpareto(n=n, d=d, model="HR", par=G)
  
  
  if(reg_method=="emst"){
    ptm <- proc.time()[1]
    fit <- emst(data = X, p=p, method = "vario")
    time <- proc.time()[1] - ptm
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=fit$graph))
  }
  else if(reg_method=="MTP2"){
    G_emp <- emp_vario(data = X, p = p)
    ptm <- proc.time()[1]
    G_emtp2 <- emtp2(G_emp, tol=1e-6,verbose = FALSE)$G_emtp2 
    time <- proc.time()[1] - ptm
    adj_emtp2 <- (abs(Gamma2Theta(G_emtp2)) >= 1e-4) 
    graph_emtp2 <- igraph::graph_from_adjacency_matrix(adj_emtp2, mode = "undirected", diag = FALSE)
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=graph_emtp2))
  }
  else{
    ptm <- proc.time()[1]
    fit <- eglearn2(data = X, p=p, rholist = rholist, reg_method = reg_method, ic = ic)
    time <- proc.time()[1] - ptm
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=fit$graph[[i]]))
    if(reg_method=="ns" | ic=="hr") F1_ic <- sapply(1:3, FUN = function(i) F1_score(g=g, gest=fit$graph_ic[[i]])) # 
  }
  
  if(reg_method=="ns" | reg_method=="glasso"){
    if(any(unlist(lapply(fit$graph, is_connected)))) F1_max <- max(F1[unlist(lapply(fit$graph, is_connected))])
    else F1_max <- 0
  }
  else F1_max <- max(F1)
  
  tbl <- tibble(type = paste0("time"), 
                value = time) %>% 
    bind_rows(tibble(type = paste0("F1_score", 1:length(rholist)),
                     value =  F1)) %>% 
     bind_rows(tibble(type = c("incoh_pos_GLi", "incoh_pos_NSi"),
                      value =  incoh_pos)) %>% 
    bind_rows(tibble(type = c("partial_pos", "partial_neg"),
                     value =  c(partial_pos,partial_neg))) %>% 
     bind_rows(tibble(type = c("F1_max"),
                      value =  F1_max)) %>% 
    bind_rows(tibble(type = c("F1_ns_aic", "F1_ns_bic", "F1_ns_mbic"),
                     value =  F1_ic))
  
  return(tbl)
}



F1_score <- function(g, gest) {
  (2*ecount(intersection(gest, g)))/( 2*ecount(intersection(gest, g)) + 
                                        ecount(intersection(complementer(g), gest)) +
                                        ecount(intersection(g, complementer(gest))))
}



incoherence_pess <- function(Sig, S, Sc){
  if(length(Sc)!=0)
    1 - max(sapply(Sc, function(e) sum(abs( Sig[e, S] %*% solve(Sig[S, S]) ))))
  else
    1
}
  
incoherence <- function(Sig, S, Sc, signs){
  if(length(Sc)!=0)
    1 - max(abs(sapply(Sc, function(e) sum( Sig[e, S] %*% solve(Sig[S, S]) %*% signs) )))
  else
    1
}
  

GLNSi <- function(Sig, cor.scale=TRUE, pessimistic=FALSE){
  
  d <- ncol(Sig)
  if(cor.scale) Sig <- ( diag(Sig)^(-0.5) * diag(d) ) %*% Sig %*% ( diag(Sig)^(-0.5) * diag(d) )
  
  Th <- solve(Sig)
  Th <- Th * (abs(Th) >= 1e-6)
  
  gr <- graph_from_adjacency_matrix(Th != 0, mode = "undirected")
  S <- which(as.vector(gr[1:d]) == 1)
  Sc <- which(as.vector(gr[1:d]) == 0)
  H <- kronecker(Sig, Sig)
  GLi <- incoherence_pess(H, S, Sc)
  
  NSi <- min(sapply(1:d, function(a){
    if(sum(Th[-a,a] != 0) == 0) return(1)
    else{
      S <- which(Th[,a] != 0)
      S <- S[S != a]
      Sc <- which(Th[,a] == 0)
      if(pessimistic) return(incoherence_pess(Sig, S, Sc))
      else{
        signs <- sign(Th[S,a])
        return(incoherence(Sig, S, Sc, signs))
      }
    }
  }))
  #GLi=NA
  return(c(GLi, NSi))
  
}

# If pessimistic = F (the default), NSi is equal to the pessimistic lower bound on the incoherence,
# with the operator norm
# Otherwise NSi is the exact incoherence depending on the signs of the regression coefficients
# GLi is unchanged by pessimistic
Gamma2Inc <- function(G, cor.scale=TRUE, pessimistic=FALSE){
  
  
  Inc <- t(sapply(1:nrow(G), function(m){
    
    Sig <- Gamma2Sigma(G, m)
    return(GLNSi(Sig, cor.scale, pessimistic))
    
  }))
  
  return(list(GLi=Inc[,1], NSi=Inc[,2]))
  
}


# traditional criteria
aic <- function(n, p) 2
bic <- function(n, p) log(n)

# modified BIC of Wang & Leng, JRSSB 2009
mbic <- function(n, p) log(n) * log(log(p))



wrapper_sim <- function(i, rowid, sim_fn, sim_fn_args){
  ## apply arguments sim_fn_args[i] to sim_fn
  ## Ahttps://uniart1.wixsite.com/uni-artrgs:
  ##    - i (integer): row to consider from sim_fn_args
  ##    - rowid (integer): unique identifier of the current simulation row
  ##      (not necessarily equal to i)
  ##    - sim_fn (function): function to run
  ##    - sim_fn_args (tibble): tibble with arguments to pass to sim_fn
  ##
  ## Returns:
  ##    - tibble with simulation results
  
  do.call(what = sim_fn, args = sim_fn_args[i, ]) %>% #ML: fixed the name of the last argument. It used to be fun_args
    mutate(rowid = rowid)
}

set_rng <- function(tbl, seed){
  ## adds to tbl a column with seeds to generate independent streams of random
  ## numbers.
  ##
  ## Args:
  ##     - tbl: a tibble where the columns contain the parameter settings and the
  ##       rows contain the simulation runs.
  ##
  ## Returns:
  ##     The function returns tbl appending a column with seeds used to generate
  ##     independent streams of random numbers.
  ##
  ## Note:
  ##     This function creates ensures that the simulations are fully repeatable.
  ##     This is possible because it assigns to each simulation run a unique 
  ##     random seed (generated with L'Ecuyer RNG method, which is suitable 
  ##     for parallel processing, too).
  
  m <- n_groups(tbl)
  group_idxs <- group_indices(tbl)
  
  # create independent RNG streams with L'Ecuyer method
  rng <- RNGseq(m, seed = seed, simplify = FALSE)
  rng <- rng[group_idxs]
  
  # add RNG streams to tbl
  tbl$rng <- rng
  
  # return tibble
  return(tbl)
}

rep_tibble <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
} 

assign_random_seed <- function(tbl, grouping_vars, seed){
  ## tibble character_vector integer -> tibble
  ## assign random seed according to the variables in grouping_vars
  if (is.null(grouping_vars)){
    tbl <- tbl %>%
      rowwise()
  } else {
    tbl <- tbl %>% 
      group_by(across(all_of(grouping_vars)))
  }
  
  tbl %>% 
    set_rng(seed) %>% 
    ungroup()
}


#rep_tibble_new solves an issue with package intersection

rep_tibble_new <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
}



generate_BA_model <- function(d,m){
  g <- sample_pa(n=d, m=m, zero.appeal=1,directed=F)
  W <- as_adj(g, sparse=F) * matrix(runif(d^2,2, 5), nrow=d)
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  O <- diag(rowSums(W)) - W
  G <- Theta2Gamma(O)
  return(list(G = G, graph = Gamma2graph(Gamma=G)))
}



generate_block_model <- function(ncliques, clique_size, alphad = 1){
  kk <- clique_size
  GG <- matrix(NA, ncliques*(kk-1) + 1, ncliques*(kk-1) + 1)
  for(i in 1:ncliques){
    bigS <- rcorrmatrix(kk, alphad = alphad)
    G1 <- Sigma2Gamma(bigS, full = TRUE)
    if(i==1) GG[1:kk, 1:kk] <- G1
    else GG[(kk + (i-2)*(kk-1)):(kk + (i-2)*(kk-1) + kk - 1), (kk + (i-2)*(kk-1)):(kk + (i-2)*(kk-1) + kk - 1)] <- G1
  }
  G <- complete_Gamma(GG)
  round(Gamma2Theta(G),2)
  sum(sum(round(Gamma2Theta(G),2) > 0)) - nrow(G)
  sum(sum(round(Gamma2Theta(G),2) < 0))
  
  return(list(G=G, graph = Gamma2graph(G)))
}

save_myplot <- function(plt, plt_nm,
                        width, height, 
                        width_pdf = 50, height_pdf = 50,
                        crop = TRUE, cairo = TRUE) {
  
  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }
  
  if (cairo) {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"), 
           device = cairo_pdf, family = "Arial")
  } else {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"))
  }
  
  if (crop){
    knitr::plot_crop(plt_nm)
  } 
}



lst_methods <- list("emst" = "emst",
                    "MTP2" = "emtp2",
                    "glasso" = "gl",
                    "ns" = "ns",
                    "F1_ns_aic" = "ns_aic",
                    "F1_ns_bic" = "ns_bic",
                    "F1_ns_mbic" = "ns_mbic"
)

my_palette <- list(
  "red" = "#D55E00",
  "blue" = "#0072B2", 
  "green" = "#009E73",
  "yellow" = "#E69F00",
  "pink" = "#CC79A7",
  "light_blue" = "#56B4E9",
  "grey" = "#999999",
  "background" = "#332288"
)



# my_palette_methods <- list(
#   c("reg_method" = lst_methods$emst, "color" = my_palette$red, "fill" = "white"),
#   c("reg_method" = lst_methods$glasso, "color" = my_palette$blue, "fill" = "white"),
#   c("reg_method" = lst_methods$MTP2, "color" = my_palette$green, "fill" = "white"),
#   c("reg_method" = lst_methods$ns, "color" = my_palette$yellow, "fill" = "white"),
#   c("reg_method" = lst_methods$F1_ns_aic, "color" = "black", "fill" = my_palette$pink),
#   c("reg_method" = lst_methods$F1_ns_bic, "color" = "black", "fill" = my_palette$light_blue),
#   c("reg_method" = lst_methods$F1_ns_mbic, "color" = "black", "fill" = my_palette$grey)
# ) %>% 
#   purrr::transpose() %>% 
#   as_tibble() %>% 
#   unnest(cols = c(reg_method, color, fill)) 

# my_col <-  my_palette_methods %>% 
#   dplyr::select(reg_method, color) %>% 
#   deframe()

# my_fill <-  my_palette_methods %>% 
#   dplyr::select(reg_method, fill) %>% 
#   deframe()


refactor_methods <- function(methods, lst_methods){
  ## character_vector list with mapping -> factor
  ## refactor column with methods
  
  lst_methods <- lst_methods
  
  
  unique_methods <- unique(methods)
  
  new_levels <- names(lst_methods)
  new_labels <- lst_methods %>% unlist() %>% unname()
  
  factor(methods,
         levels = new_levels,
         labels = new_labels)
}


theme_fct <- function(font_size1=11,  font_size2=7.5){
  theme_set(theme_bw() +
              theme(
                plot.background = element_blank(),
                panel.background = element_blank(),
                legend.background = element_blank(),
                strip.background = element_rect(fill = "white"),
                plot.caption=element_text(size=font_size2, hjust=0, 
                                          margin=margin(t=15)),
                text = element_text(size = font_size1),
                axis.ticks = element_blank(),
                axis.text = element_text(size = font_size1),
                panel.grid.major = element_line(size = 0.25)
              ) 
            #+
            # theme_cowplot(font_size = 11)
  )
}

theme_fct()


create_palette_levels <- function(reg_method_n, palette_tbl){
  
  str_spl <- strsplit(reg_method_n, "__")
  
  my_tbl <- tibble(
    reg_method = purrr::map_chr(str_spl, function(el){el[1]}),
    level =  purrr::map_chr(str_spl, function(el){el[2]})
  ) %>%
    left_join(palette_tbl, by = "reg_method") %>%
    mutate(fill = if_else(level == "2", color, fill)) %>%
    mutate(reg_method_lev = paste(reg_method, level, sep = "__"))
  
  my_col <-  my_tbl %>%
    dplyr::select(reg_method_lev, color) %>%
    deframe()
  
  my_fill <-  my_tbl %>%
    dplyr::select(reg_method_lev, fill) %>%
    deframe()
  
  
  list(cols = my_col, fills = my_fill)
  
}
