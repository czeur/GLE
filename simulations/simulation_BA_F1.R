# library(graphicalExtremes)
# library(igraph)
library(parallel)

source("functions/sim_graph.R")

# Path to store figures
figpath <- "figures/"
respath <- "results/"
num_cores <- 10
set.seed(100)

# For each d value, there is corresponding
# n: sample size
# q.threshold: choice of threshold
# nsim: number of replication

dset <- c(100)
qset <- c(1, 2)

nsim <- 100

for (d in dset){
    if (d == 20){
        lambda.range <- 10 ^ seq(from = -2, to = 0, by = 0.1)
    } else if (d == 100){
        lambda.range <- 10 ^ seq(from = -2.4, to = -1.4, by = 0.2)
    }
    for (q in qset){
        if (q == 1){
            kset <- c(1, 2.5, 5) * d
        } else if (q == 2){
            kset <- c(1, 5, 10) * d
        }
        graph <- sample_pa(n = d, m = q, zero.appeal = 1, directed = FALSE)
        graph.true <- graph
        W_mat <- as_adj(graph, sparse = FALSE) * matrix(runif(d^2, 2, 5), nrow = d)
        W_mat[lower.tri(W_mat)] <- t(W_mat)[lower.tri(W_mat)]
        Theta <- diag(rowSums(W_mat)) - W_mat
        Gamma <- Theta2Gamma(Theta)
        cat("\n Finish setting up Theta for: d = ", d, " q = ", q, "\n")

        for (k in kset){
            n <- ceiling(k ^ (1 / 0.7))
            q.threshold <- 1 - k / n

            all_objects <- ls() # Use one less core than available
            cl <- makeCluster(num_cores)
            clusterExport(cl, all_objects)
            clusterEvalQ(cl, library(graphicalExtremes))
            clusterEvalQ(cl, library(glmnet))
            clusterEvalQ(cl, library(igraph))

            cat("\n Ready to start simulation with n = ", n, " k = ", k, "\n")
            start_time <- Sys.time()
            
            total_jobs <- nsim
            batch_size <- num_cores
            
            for (i in seq(1, total_jobs, by = batch_size)) {
                print(i)
                jobs_to_run <- min(batch_size, total_jobs - i + 1)
                random_seeds <- sample.int(1e6, jobs_to_run)
                clusterExport(cl, "random_seeds")

                # results <- c(results, mclapply(1:jobs_to_run, function(j){
                #     gc()
                #     library(graphicalExtremes)
                #     library(glmnet)
                #     library(igraph)
                #     set.seed(random_seeds[j])
                #     sim_graph(Theta, Gamma, n, d, lambda.range = lambda.range, q.threshold, thres = 1e-3, iter.max = 1e3, BIC_method = "maxMBIC")
                # }, mc.cores = 8))
                
                
                parLapply(cl, 1:jobs_to_run, function(j) {
                    set.seed(random_seeds[j])
                    fit <- sim_graph(Theta, Gamma, n, d, lambda.range = lambda.range, q.threshold, thres = 1e-2, iter.max = 1e3, BIC_method = "MBIC")
                     fn =  paste(respath,"BA_q" , q, "/Results_d_", d, "_n_", n, "_k_", k, "_", random_seeds[j], ".Rdata", sep = "")
                    save(fit, Theta, file = fn)
                    rm(fit)
                    gc()
                })
                gc()
            }
            stopCluster(cl)
            cat("\n Time used for simulation: ", Sys.time() - start_time)
            cat("\n")
            # fn =  paste(respath,"BA_q" , q, "/Results_d_", d, "_n_", n, "_k_", k, ".Rdata", sep = "")
            # save(results, Theta, file = fn)
        }
    }
}