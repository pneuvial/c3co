## a simple problem with very little noise
nbBkps <- 10
nbClones <- 4

Z <- matrix(1, nrow = nbClones, ncol=nbBkps+1)
Z[2, 2] <- 2
Z[3, 5:6] <- 2
Z[4, 9:10] <- 2
matplot(t(Z), t='s')

W <- diag(rep(1, nrow(Z)))
W <- rbind(W, W)
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0.05), 
            nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E
round(Y, 2)
matplot(t(Y), t='s')  ## pretty low noise

lambda <- 1e-4
params <- list(lambda=lambda, nb.arch=2:nrow(W))
res <- fitC3co(Y, parameters.grid=params)
fails <- sapply(res$fit, function(x) x@failure)
print(fails)

res$config$best
pvePlot2(res$config$best)
plot(res$config$best$BIC)
plot(res$config$best$logLik)
diff(res$config$best$logLik)
pen <- res$config$best$logLik - res$config$best$BIC
plot(pen)

## with another lambda ? (testing sensitivity of PVE)
lambdas <- 10^(-(10:4))
stats <- lapply(lambdas, FUN=function(lambda) {
    params <- list(lambda=lambda, nb.arch=2:nrow(W))
    res <- fitC3co(Y, parameters.grid=params)
    fails <- sapply(res$fit, function(x) x@failure)
    if (any(fails)) {
        print(lambda)
        print(fails)
    }
    res$config$best
})
stats <- Reduce(rbind, stats)

library("ggplot2")
ggplot(stats, aes(x=nb.feat, y=PVE, color=lambda1, group=lambda1)) + 
    geom_line()
ggplot(stats, aes(x=nb.feat, y=BIC, color=lambda1, group=lambda1)) + 
    geom_line() 
ggplot(stats, aes(x=nb.feat, y=logLik, color=lambda1, group=lambda1)) + 
    geom_line() 
ggplot(stats, aes(x=nb.feat, y=log(loss), color=lambda1, group=lambda1)) + 
    geom_line() 

