source('glasso_c_function.R')

d <- 100
A <- matrix(rnorm(d^2),d,d)
S <- A%*%t(A) + t(A)%*%A

lambda <- 10
c <- 1/(d*eigen(S)$values[1])

test <- glasso_c(S,lambda=lambda,c=c)
Theta <- test$Theta-c
View(Theta)
