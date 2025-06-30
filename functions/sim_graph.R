## Using BIC to tune

source("functions/BIC.R")

# Parameters
# Theta: the true Theta matrix used in the simulation
# Gamma: the true Gamma matrix used in the simulation
# n: sample size
# d: dimension

# lambda.range, a range of lambda to be chosen
# Setting lambda.range = 0 will select automatically a range

# q.threshold: choice of threshold used in estimation
# thres: threshold used to trunctate to zero after glasso_c
# iter.max: max iteration used in glasso_c

# BIC_method: BIC used in the calculaton, options are
# "mean": average BIC from d graphs
# "max": maximum BIC from d graphs
# "MBIC": total MBIC
# "maxMBIC": maximum MBIC

sim_graph <- function (Theta, Gamma, n, d, lambda.range = 0, q.threshold, thres, iter.max, BIC_method = "maxMBIC"){

    data <- rmpareto(n=n, model="HR", d=d, par=Gamma)

    fit <- glasso_bic(data, lambda.range = lambda.range , q.threshold = q.threshold,
    thres = thres, iter.max = iter.max, BIC_method = BIC_method)

    # rm(list = "data")
    
    # gc(verbose = FALSE)
    
    return(fit)
}


sim_graphs_multilambda = function(Theta, n, d, q.threshold, lambda.range = 0, thres, iter.max, BIC_method = "maxMBIC"){

    data <- rmpareto(n=n, model="HR", d=d, par=Gamma) 

    W <- est_W(data, q.threshold)
    Gamma <- emp_vario(data, p = q.threshold)
    c <- 1 / (d * eigen(W$cov)$values[1])

    k <- n * (1 - q.threshold)

    results <- list()

    nlambda <- length(lambda.range)
    connected <- rep(FALSE, nlambda)
    IC <- rep(NA, nlambda)

    for (i in 1: nlambda){
      lambda <- lambda.range[i]

      res <- glasso_c_reest(W$cov, lambda, c, thres, iter.max)
      results <- c(results, list(res$Theta))
      
      graph <- igraph::graph_from_adjacency_matrix(!(res$graph), diag = FALSE,mode = "undirected")
      if(!is.connected(graph)){
          next
      } else {
            connected[i] <- TRUE
            Gamma_hat <- complete_Gamma(Gamma = Gamma, graph = graph)
            Theta <- Gamma2Theta(Gamma_hat)

            BIC <- rep(NA, d)
            for (j in 1:d){
                Theta_part <- Theta[-j, -j]
                missing <- sum(res$graph[-j,-j])
                p <- ((d - 1) * d / 2 - missing / 2 )
                
                if (BIC_method == "mean" || BIC_method == "max"){
                    BIC[j] <- - log(det(Theta_part)) + sum(diag(W$subcovlist[[j]] %*% Theta_part)) + p * log(k) / k
                    # BIC[j] <- - alldet + sum(diag(W$subcovlist[[j]] %*% Theta_part)) + p * log(k) / k
                } else if (BIC_method == "MBIC" || BIC_method == "maxMBIC" || BIC_method == "medianMBIC") {
                BIC[j] <- - log(det(Theta_part)) + sum(diag(W$subcovlist[[j]] %*% Theta_part)) + p * log(k) * log(log(d - 1)) / k
                    # BIC[j] <- - alldet + sum(diag(W$subcovlist[[j]] %*% Theta_part)) + p * log(k) * log(log(d - 1)) / k
                }
            }
            if (BIC_method == "max" || BIC_method == "maxMBIC"){
                IC[i] <- max(BIC)
            } else if (BIC_method == "mean" || BIC_method == "MBIC") {
                IC[i] <- mean(BIC)
            } else if (BIC_method == "medianMBIC"){
                IC[i] <- quantile(BIC,0.5)
            }
        }
    }
    
    return(list(results, connected, IC))
}