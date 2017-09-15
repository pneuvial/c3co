library(dplyr)
## Example 1 with noise
## a simple problem with noise
nbBkps <- 10
nbClones <- 4

Z <- matrix(1, nrow = nbClones, ncol=nbBkps+1)
Z[2, 2] <- 2
Z[3, 5:6] <- 2
Z[4, 9:10] <- 2
matplot(t(Z), t='s')
W <- diag(rep(1, nrow(Z)))
W <- rbind(W, W)
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0.1), 
            nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E
round(Y, 2)
matplot(t(Y), t='s')  ## pretty low noise

lambda <- 1e-4
params <- list(lambda=lambda, nb.arch=2:nrow(W))
resPosFus <- c3co::positiveFusedLasso(list(Y), list(t(Z)), lambda) 
res <- fitC3co(Y, parameters.grid=params)
res$config$best
pvePlot2(res$config$best)
plot(res$config$best$BIC)
plot(res$config$best$logLik)
diff(res$config$best$logLik)
pen <- res$config$best$logLik - res$config$best$BIC
plot(pen)
res$fit[[3]]@S$Z %>% round(2) %>% t
par(mfrow=c(2,1), mar=c(2,5,2,5))
matplot(res$fit[[3]]@S$Z %>% round(2), type="s", ylab="Z")
matplot(t(Z), type="s", ylab="Z true")
par(mfrow=c(2,1))
matplot(res$fit[[2]]@E$Y %>% t, type="s")
matplot(t(Y), type="s")
res$fit[[3]]@W %>% heatmap.3()
## - - - - - - - - - - - - - - - - - - - - - -
## noisy example without normal
## - - - - - - - - - - - - - - - - - - - - - -
lambda <- 1e-4
W <- diag(rep(1, nrow(Z)))
W <- W[-1, ]
W <- rbind(W, W)
stopifnot(all(rowSums(W)==1))
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0.1), 
            nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E
lambda <- 1e-4
params <- list(lambda=lambda, nb.arch=2:nrow(W))
res <- fitC3co(Y, parameters.grid=params)
res$config$best

pvePlot2(res$config$best)
par(mfrow=c(2,1), mar=c(2,5,2,5))
matplot(res$fit[[2]]@S$Z %>% round(2), type="s", ylab="Z")
matplot(t(Z), type="s", ylab="Z true")
par(mfrow=c(2,1))
matplot(res$fit[[2]]@E$Y %>% t, type="s")
matplot(t(Y), type="s")
res$fit[[2]]@W %>% heatmap.3()

## - - - - - - - - - - - - - - - - - - - - - - - - - - -
## small example normal contamination but normal clone *not observed*
## - - - - - - - - - - - - - - - - - - - - - - - - - - -
W <- diag(rep(1, nrow(Z)))
lambda <- 1e-4
W1 <- W[-1, ]
stopifnot(all(rowSums(W1)==1))
E <- matrix(rnorm(nrow(W1)*ncol(Z), sd=0.1), nrow=nrow(W1), ncol=ncol(Z))
Y <- W1 %*% Z + E
lambda <- 1e-3
params <- list(lambda=lambda, nb.arch=2:3)
res <- fitC3co(Y, parameters.grid=params)

res$config$best
pvePlot2(res$config$best)
par(mfrow=c(2,1), mar=c(2,5,2,5))
matplot(res$fit[[2]]@S$Z %>% round(2), type="s", ylab="Z")
matplot(t(Z), type="s", ylab="Z true")
par(mfrow=c(2,1))
matplot(res$fit[[2]]@E$Y %>% t, type="s")
matplot(t(Y), type="s")
res$fit[[2]]@W %>% heatmap.3()

## - - - - - - - - - - - - - - - - - - - - - -
## non-totally trivial W and normal contamination 
## Noised observations
## Define W by rSparseWeightMatrix function in c3co package
## - - - - - - - - - - - - - - - - - - - - - -
W <- rSparseWeightMatrix(nb.samp=7, nb.arch=3, sparse.coeff=0.5)
W <- (1-rowSums(W %>% as.matrix)) %>% cbind(W %>% as.matrix)
stopifnot(all(rowSums(W)==1))
Z <- matrix(1, nrow = nbClones, ncol=nbBkps+1)
Z[2, 2] <- 2
Z[3, 5:6] <- 2
Z[4, 9:10] <- 2
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0.1), nrow=nrow(W), ncol=ncol(Z))
Y <- W  %*% Z + E
matplot(t(Y), type="s")
lambda <- 1e-4
params <- list(lambda=lambda, nb.arch=2:6)
res <- fitC3co(Y, parameters.grid=params )
pvePlot2(res$config$best)
res$fit[[3]]@S$Z1 %>% t %>% round(1)
res$fit[[3]]@W %>% round(2)
W %>% round(2) %>% heatmap.3(cellnote=W)
res$fit[[3]]@W %>% round(2) %>% heatmap.3(cellnote=res$fit[[3]]@W  %>% round(2) )

res$fit[[3]]@S$Z %>% round(2) %>% t
par(mfrow=c(2,1), mar=c(2,5,2,5))
matplot(res$fit[[3]]@S$Z %>% round(2), type="s", ylab="Z")
matplot(t(Z), type="s", ylab="Z true")
par(mfrow=c(2,1))
matplot(res$fit[[3]]@E$Y %>% t, type="s")
matplot(t(Y), type="s")
