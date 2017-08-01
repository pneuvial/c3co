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
lambda <- 1e-5
W <- diag(rep(1, nrow(Z)))
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E

res <- positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(res@W, 2)  ## W
round(t(res@S$Z), 2)  ## Z
round(res@mu, 2) ## 0 because at initialization, getW gets the right W straight ahead.

## - - - - - - - - - - - - - - - - - - - - 
## with the correct Z as a starting point
## but without the normal clone !
## - - - - - - - - - - - - - - - - - - - - 
res <- positiveFusedLasso(list(Y), list(t(Z[-1, ])), lambda) ## with the correct Z minus normal subclone
round(res@W, 2)       ## W[-1, ] OK but W[1, ] not OK !!
round(t(res@S$Z), 2)  ## Z is okay but already not perfect
round(res@E$Y, 2) 
round(res@mu, 2)      ## mu[1]!=0 because it tries to reconstruct Y[1, ] (normal patient)


## what we would like in this situation is W[1, ] equal to 0, but this is not possible due to the convexity constraint on W! => something is wrong

## a sanity check
Y2 <- res@W %*% t(res@S$Z)  ## without adding back mu
Y3 <- sweep(res@E$Y, 1, res@mu, "-") ## removing mu again (equivalent)
identical(Y2, Y3)  ## TRUE

## - - - - - - - - - - - - - - - - - - - - 
## what if normal clone not in Z0 
## and not in Y either? 
## Back to perfection!
## - - - - - - - - - - - - - - - - - - - - 
res <- positiveFusedLasso(list(Y[-1, ]), list(t(Z[-1, ])), lambda) ## with the correct Z,Y= minus normal subclone,sample
round(res@W, 2)
round(t(res@S$Z), 2)
round(res@E$Y, 2) 
round(res@mu, 2) 

## so the pb is to reconstruct a normal y w/o an explicit normal subclone 

## - - - - - - - - - - - - - - - - - - - - 
## non-totally trivial W
## - - - - - - - - - - - - - - - - - - - - 
lambda <- 1e-4
W <- rbind(c(1, 1, 0, 0)/2,
           c(1, 0, 1, 0)/2,
           c(1, 0, 0, 1)/2)
stopifnot(all(rowSums(W)==1))
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E

res <- positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(res@W, 2) - W   ## good
round(t(res@S$Z), 1)  ## not good! 
round(t(res@S$Z), 1) - Z
round(res@mu, 2)
Y-res@E$Y  ## still off by a constant!

modelFitStats(res)  ## OMG we're back to negative PVE :(

