## a simple problem with very little noise
nbBkps <- 10
nbClones <- 4

Z <- matrix(1, nrow = nbClones, ncol=nbBkps+1)
Z[2, 2] <- 2
Z[3, 5:6] <- 2
Z[4, 9:10] <- 2
matplot(t(Z), t='s')

## - - - - - - - - - - - - - - - - - - - - 
## with the correct Z as a starting point
## - - - - - - - - - - - - - - - - - - - - 
lambda <- 1e-3
W <- diag(rep(1, nrow(Z)))
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E

res <- c3co::positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(res@W, 2)  ## W
round(t(res@Zt$Z), 2)  ## Z

## - - - - - - - - - - - - - - - - - - - - 
## with the correct Z as a starting point
## but without the normal clone !
## - - - - - - - - - - - - - - - - - - - - 
res <- c3co::positiveFusedLasso(list(Y), list(t(Z[-1, ])), lambda) ## with the correct Z minus normal subclone
round(res@W, 2)       ## W[-1, ] OK but W[1, ] not OK !!
round(t(res@Zt$Z), 2)  ## Z is almost not perfect
round(res@E$Y, 2) 

## what we would like in this situation is W[1, ] equal to 0, but this is not possible due to the convexity constraint on W! Conclusion: force a column of 1 in Z as input of positiveFusedLasso

## - - - - - - - - - - - - - - - - - - - - 
## what if normal clone not in Z0 
## and not in Y either? 
## Back to perfection!
## - - - - - - - - - - - - - - - - - - - - 
res <- c3co::positiveFusedLasso(list(Y[-1, ]), list(t(Z[-1, ])), lambda) ## with the correct Z,Y= minus normal subclone,sample
round(res@W, 2)
round(t(res@Zt$Z), 2)
round(res@E$Y, 2) 

## so the pb is to reconstruct a normal y w/o an explicit normal subclone 

## - - - - - - - - - - - - - - - - - - - - 
## what if normal clone in Z0 
## and not in Y either? 
## Manage this internally by removing the column of Z
## This can be detected by checking the rank in W so that WtW remain invertible
## - - - - - - - - - - - - - - - - - - - - 
res <- c3co::positiveFusedLasso(list(Y[-1, ]), list(t(Z)), lambda) ## with the correct Z,Y= minus normal subclone,sample
round(res@W, 2)
round(t(res@Zt$Z), 2)
round(res@E$Y, 2) 

## - - - - - - - - - - - - - - - - - - - - - -
## non-totally trivial W: normal contamination
## - - - - - - - - - - - - - - - - - - - - - -
lambda <- 1e-4
W <- rbind(c(1, 0, 0, 0),
           c(1, 1, 0, 0)/2,
           c(1, 0, 1, 0)/2,
           c(1, 0, 0, 1)/2)
stopifnot(all(rowSums(W)==1))
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E

res <- c3co::positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(res@W, 2)
round(t(res@Zt$Z), 1)
Y-res@E$Y  ## rather good

modelFitStats(res)  ## no more negative PVE :>)

## - - - - - - - - - - - - - - - - - - - - - - - - - - -
## normal contamination but normal clone *not observed*
## - - - - - - - - - - - - - - - - - - - - - - - - - - -
lambda <- 1e-4
W1 <- W[-1, ]
stopifnot(all(rowSums(W1)==1))
E <- matrix(rnorm(nrow(W1)*ncol(Z), sd=0), nrow=nrow(W1), ncol=ncol(Z))
Y <- W1 %*% Z + E

res <- c3co::positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(res@W, 2)
round(t(res@Zt$Z), 1)
Y-res@E$Y  ## rather good given that we have 3 patients and 4 archetypes !!!
## in fact here the *true* W is singular because nrow(W) > ncol(W) !

modelFitStats(res)  ## no more negative PVE :>)


## - - - - - - - - - - - - - - - - - - - - - - - - - - -
## normal contamination 
## but normal clone *not observed*
## but nb of patients >= nb of archetypes (otherwise can't expect recovering truth)
## - - - - - - - - - - - - - - - - - - - - - - - - - - -
lambda <- 5e-4
W <- rbind(c(1, 1, 0, 0)/2,
           c(1, 0, 1, 0)/2,
           c(1, 0, 0, 1)/2,
           c(0, 0, 1, 1)/2)
qr(W)$rank  ## not singular
stopifnot(all(rowSums(W)==1))
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E

res <- c3co::positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(res@W, 2)## First line is not perfect, due to subclone 2 only in profile 1
round(t(res@Zt$Z), 2)## Subclones are quite good except subclone 2
Y-res@E$Y  ## rather good
modelFitStats(res)  


## - - - - - - - - - - - - - - - - - - - - 
## Colinearity test
## - - - - - - - - - - - - - - - - - - - - 
nbBkps <- 10
nbClones <- 4
Z <- matrix(1, nrow = nbClones, ncol=nbBkps+1)
Z[2, 2] <- 2
Z[3, 5:6] <- 2
Z[4,] <- 0.5*(Z[2,]+Z[3,])
matplot(t(Z), t='s')
## - - - - - - - - - - - - - - - - - - - - 
## with the correct Z as a starting point and no noise
## - - - - - - - - - - - - - - - - - - - - 
lambda <- 1e-3
W <- rbind(c(0, 1, 0, 0),
           c(0, 1, 1, 0)/2,
           c(0, 1, 0, 1)/2,
           c(0, 0, 1, 1)/2)
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0.0), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E
res <- c3co::positiveFusedLasso(list(Y), list(t(Z[-c(1), ])), lambda) 
round(res@W, 2)  ## W
round(t(res@Zt$Z), 2)  ## Z
qr(t(res@Zt$Z))$rank 
par(mfrow=c(2,2))
matplot(round(res@Zt$Z, 2), type="s")
matplot(t(Z[-1,]), type="s")
matplot(t(Y), type="s")
matplot(t(res@E$Y), type="s")
heatmap.3(round(res@W, 2) )
heatmap.3(W[,-1]) 

## - - - - - - - - - - - - - - - - - - - - 
## with the correct Z as a starting point with noise
## - - - - - - - - - - - - - - - - - - - - 
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0.1), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E
res <- c3co::positiveFusedLasso(list(Y), list(t(Z[-c(1), ])), lambda=1e-4)
round(res@W, 2)  ## W
round(t(res@Zt$Z), 2)  ## Z
qr(t(res@Zt$Z))$rank 
par(mfrow=c(2,2))
matplot(round(res@Zt$Z, 2), type="s")
matplot(t(Z[-1,]), type="s")
matplot(t(Y), type="s")
matplot(t(res@E$Y), type="s")
heatmap.3(round(res@W, 2) )
heatmap.3(W[,-1]) 
