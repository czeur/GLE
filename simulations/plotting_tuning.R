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


d <- 20
qset <- c(1, 2)

nsim <- 100

lambda.range <- seq(from = -2, to = 0, by = 0.1)

nlambda <- length(lambda.range)

for (q in qset){
    
    cat("\n Start for: d = ", d, " q = ", q, "\n")

    F1 <- matrix(0, nsim, nlambda)
    IC <- F1
    connected <- F1

    isopt <- matrix(FALSE, nsim, nlambda)

    
    fn =  paste("results/Multiple tuning_d_", d, "_q_", q, ".Rdata", sep = "")
    load(file = fn)

    for (i in 1:nlambda){
        for (j in 1:nsim){
            F1[j, i] <- F1score(Theta, results[[j]][[1]][[i]])
        }
    }

    for (j in 1:nsim){
        tempIC <- results[[j]][[3]]
        tempIC[is.na(tempIC)] <- Inf
        IC[j, ] <- tempIC
        connected[j, ] <- results[[j]][[2]]

        isopt[j, which.min(tempIC)] = TRUE
    }
    
    F1 <- colMeans(F1)
    connected <- colMeans(connected)

    isopt <- colMeans(isopt)

    # IC <- colMeans(IC, na.rm = TRUE)
    fn = paste("figures/Tuning_plots_q_", q, ".pdf", sep = "")
    lambda.range <- seq(from = -2, to = 0, by = 0.1)
    pdf(fn, width = 8, height = 5)
    par(mar = c(5, 4, 4, 4))
    plot(lambda.range, F1, type = "o", xlab = expression(log[10](gamma)), ylab = "F1 Score", ylim = c(min(F1), max(F1)))

    par(new=TRUE)
    plot(isopt, type = "l", lty = 2, xlab = "", axes = FALSE, ylab="")
    axis(side=4)
    mtext("Fraction of being optimal", side = 4, line = 2)
    dev.off()
}
