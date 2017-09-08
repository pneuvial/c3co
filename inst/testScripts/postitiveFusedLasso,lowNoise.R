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

res <- c3co::positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(res@W, 2)  ## W
round(t(res@S$Z), 2)  ## Z

## - - - - - - - - - - - - - - - - - - - - 
## with the correct Z as a starting point
## but without the normal clone !
## - - - - - - - - - - - - - - - - - - - - 
res <- c3co::positiveFusedLasso(list(Y), list(t(Z[-1, ])), lambda) ## with the correct Z minus normal subclone
round(res@W, 2)       ## W[-1, ] OK but W[1, ] not OK !!
round(t(res@S$Z), 2)  ## Z is okay but already not perfect
round(res@E$Y, 2) 

## what we would like in this situation is W[1, ] equal to 0, but this is not possible due to the convexity constraint on W! => something is wrong

## - - - - - - - - - - - - - - - - - - - - 
## what if normal clone not in Z0 
## and not in Y either? 
## Back to perfection!
## - - - - - - - - - - - - - - - - - - - - 
res <- c3co::positiveFusedLasso(list(Y[-1, ]), list(t(Z[-1, ])), lambda) ## with the correct Z,Y= minus normal subclone,sample
round(res@W, 2)
round(t(res@S$Z), 2)
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
round(t(res@S$Z), 2)
round(res@E$Y, 2) 

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

res <- c3co::positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(t(res@S$Z), 1)
Y-res@E$Y  ## rather good

modelFitStats(res)  ## no more negative PVE :>)

