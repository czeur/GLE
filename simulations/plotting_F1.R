library(graphicalExtremes)
library(igraph)
source('functions/F1score.R')

# Path to store figures
figpath <- "figures/"

# Path storing results
respath <- "results/"


# Set dimensions, 200 or 20
# For each d value, there is corresponding sample size
# nsim: number of replication


dset <- c(20, 100)
qset <- c(1, 2)

nsim <- 100

for (d in dset){
    for (q in qset){
        if (q == 1){
            ratio <- c(1, 2.5, 5)
        } else if (q == 2){
            ratio <- c(1, 5, 10)
        }
        cat("\n Start for: d = ", d, " q = ", q, "\n")

        kset <- ratio * d
        F1 <- matrix(0, nsim, 3)
        lambda <- matrix(0, nsim)



        respath_q = paste(respath, "BA_q", q, sep = "")

        i <- 1
        
        for (k in kset){
            n <- ceiling(k ^ (1 / 0.7))

            file_prefix <- paste("Results_d_", d, "_n_", n, "_k_", k, sep = "")
            files_with_prefix <- list.files(path = respath_q, pattern = paste0("^", file_prefix), full.names = TRUE)

            j <- 1

            for (fn in files_with_prefix){

            # fn =  paste(respath, "Results_BA_", q, "_d_", d, "_n_", n, "_k_", k, ".Rdata", sep = "")
                load(file = fn)

                if (j == 1) Theta_check = Theta
                else {
                    if (sum(!(Theta_check == Theta)) > 0){
                        print(paste("Different Theta in", fn, "\n"))
                    }
                }

                F1[j, i] <- F1score(Theta, fit$Theta)
                lambda[j] <- fit$lambda
                j <- j + 1
            }
            i <- i + 1
        }
        x11()
        hist(lambda)
        cat("\n Minimum lambda is ", min(lambda), "; Maximum lambda is ", max(lambda), "\n")
        fn = paste("figures/F1plots_q_", q, "_d_", d, ".pdf", sep = "")
        pdf(fn, width = 8, height = 5)
        boxplot(F1, names = ratio, ylab = "F1 Scores", xlab = expression(k[n]/d))
        dev.off()
    }
}