# Function to estimate the W matrix from data

trim_matrix <- function(W, c){
    d <- dim(W)[1]
    c_max <- 1 / sum(solve(W))
    return (W - min(c, c_max) * matrix(1, d, d))
}

est_W <- function(data, q.threshold) {
    n <- dim(data)[1]
    d <- dim(data)[2]

    mdata <- matrix(0, n, d)

    for (i in 1:d){
        x <- data[, i]
        rx <- rank(x) / (n + 1)
        mdata[, i] <- 1 / (1 - rx)
    }

    cov.mat <- matrix(0, d, d) # this is S_0
    cov.k.list <- list()
    trace.sum <- 0

    for (k in 1:d){
        data.truncated <- mdata[mdata[,k]> 1 / (1 - q.threshold),]
        # data.truncated <- data[data[,k]>1,]
        w <- log(data.truncated[, -k]) - matrix(log(data.truncated[,k]), nrow(data.truncated), d-1)
        cov.k <- cov(w)
        cov.mat[-k, -k] <- cov.mat[-k, -k] + cov.k / d
        cov.k.list[[k]] <- cov.k
       trace.sum <- trace.sum + sum(cov.k) / d ^ 3
        # trace.sum <- 0
    }

    return(list(cov = trim_matrix(cov.mat, trace.sum), subcovlist = cov.k.list))
}


