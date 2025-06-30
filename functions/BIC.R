## Function for cross validation

library(caret)
library(igraph)
library(graphicalExtremes)
source("functions/glasso_c_function.R")
source("functions/est_W.R")

# Parameters
# data: input data

# lambda.range, a range of lambda to be chosen
# Setting lambda.range = 0 will select automatically a range
# That is in equal log distance grid from 1e-2 to 1e2

# q.threshold: choice of threshold used in estimation
# thres: threshold used to trunctate to zero after glasso_c
# iter.max: max iteration used in glasso_c

# BIC_method: BIC used in the calculaton, options are
# "mean": average BIC from d graphs
# "max": maximum BIC from d graphs
# "MBIC": total MBIC
# "maxMBIC": maximum MBIC

glasso_bic <- function(data, lambda.range = 0, q.threshold = 0.9,
    thres = 1e-3, iter.max = 1000, BIC_method = "maxBIC"){

    if (length(lambda.range) == 1 && lambda.range == 0){
        # lambda.range <- seq(from = 0.01, to = 0.1, by = 0.01)
        lambda.range <- seq(from = -2, to = 0, by = 0.1)
        lambda.range <- 10 ^ lambda.range
    }

    n <- dim(data)[1]
    d <- dim(data)[2]

    # Effective sample size for each graph
    k <- n * (1 - q.threshold)
    # k <- dim(data2mpareto(data, q.threshold))[1]

    nlambda <- length(lambda.range)
    IC <- rep(Inf, nlambda)

    W <- est_W(data, q.threshold)
    # Gamma <- emp_vario(data, p = q.threshold)
    # W <- Gamma2Sigma(Gamma)
    c <- 1 / (d * eigen(W$cov)$values[1])


    for (i in 1:nlambda){
        lambda <- lambda.range[i]
        res <- glasso_c_reest(W$cov, lambda, c, thres, iter.max)
        graph <- igraph::graph_from_adjacency_matrix(!(res$graph), diag = FALSE,mode = "undirected")
        if(!is.connected(graph)){
            break
        } else {

            det_Theta_part <- numeric(d)
            for (j in 1:d) det_Theta_part[j] <- det(res$Theta_hat[-j, -j])

            if (sum(det_Theta_part <= 0) > 0){
                Gamma_hat <- complete_Gamma(Gamma = Sigma2Gamma(W$cov), graph = graph)
                Theta <- Gamma2Theta(Gamma_hat)
            } else{
                Theta <- res$Theta_hat
            }
            

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
    if (min(IC) < Inf) lambda.select <- lambda.range[which.min(IC)] else lambda.select <- min(lambda.range)

    # print(cbind(lambda.range, IC))
    if (which.min(IC) == 1){
        cat("\n Warning: best lambda is at the lower bound!\n")
    } else if (which.min(IC) == nlambda){
        cat("\n Warning: best lambda is at the upper bound!\n")
    }

    # Re-estimate using the optimal lambda
    res <- glasso_c_reest(W$cov, lambda.select, c, thres, iter.max)
        # if (Theta_method == "glasso") {
        #     Theta <- res$Theta_hat
        # } else if (Theta_method == "demean") {
        #     Theta <- res$Theta_demean
        # }
    Theta <- res$Theta_hat
    Theta[res$graph] <- 0
    return(list(lambda = lambda.select, Theta = Theta))
 }
