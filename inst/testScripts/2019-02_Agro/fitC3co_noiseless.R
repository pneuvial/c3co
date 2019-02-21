len <- 100
n <- 10

z0 <- c(1, 1, 1, 1, 1)
z1 <- c(1, 0, 1, 1, 1)
z2 <- c(1, 1, 1, 0, 1)
Z <- rbind(z0, z1, z2)
nbBkps <- ncol(Z)-1

W <- diag(rep(1, nrow(Z)))  ## identity weights
Y <- W %*% Z                ## noiseless
matplot(t(Y), t='s')

params <- list(lambda=1e-5, nb.arch=2:3)
res <- fitC3co(Y, parameters.grid=params)

## 3 subclones (true value)
fit2 <- res$fit[[2]]

modelFitStats(fit2)  ## pretty good
round(fit2@W, 2)    ## W

round(fit2@Zt$Z, 2)  ## t(Z)!
round(t(fit2@Zt$Z), 2) - Z ## 

## 2 subclones (= number of non-normal subclones)
fit1 <- res$fit[[1]]
modelFitStats(fit1)  ## not so good! 

