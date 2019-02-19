## This script aims to test what it happens when we force a component to be normal.
library(dplyr)
## - - - - - - - - - - - - - - - - - - - - 
## with the correct Z as a starting point
## - - - - - - - - - - - - - - - - - - - - 
lambda <- 1e-5
Z <- matrix(1, nrow = nbClones, ncol=nbBkps+1)
Z[2, 2] <- 2
Z[3, 5:6] <- 2
Z[4, 9:10] <- 2

W <- diag(rep(1, nrow(Z)))
E <- matrix(rnorm(nrow(W)*ncol(Z), sd=0.05), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z + E

res <- positiveFusedLasso(list(Y), list(t(Z)), lambda) 
round(res@W, 2)  ## W
round(t(res@S$Z), 2)  ## Z

## should mu be 0 here??
round(res@mu, 2) ## 0
## yes because at initialization, getW gets the right Z straight ahead.
## Good results


## - - - - - - - - - - - - - - - - - - - - 
## what if normal clone not in Z0?
## - - - - - - - - - - - - - - - - - - - - 
Z0t <- initializeZt(Y, K=nrow(Z), forceNormal=FALSE)
res <- positiveFusedLasso(list(Y), Z0t, lambda) ## with the correct Z minus normal subclone
round(res@W, 2) ## shit already happens
round(t(res@S$Z), 2)  ## Z is OK (but this is a noiseless setting !!! Should be perfect)
round(res@E$Y, 2)
round(res@mu, 2) ## not 0 where the wild things are
modelFitStats(res)

Z0Nt <- initializeZt(Y, K=nrow(Z), forceNormal=TRUE)
resForceNorm <- positiveFusedLasso(list(Y), Z0Nt, lambda) 
round(resForceNorm@W, 2)
round(t(resForceNorm@S$Z), 2)  ## Z is OK (but this is a noiseless setting !!! Should be perfect)
round(resForceNorm@E$Y, 2)
round(resForceNorm@mu, 2) ## not 0 where the wild things are
modelFitStats(resForceNorm)
### Results are similar betwwen forceNormal and not forceNormal

## Test what happen with the various types of initialization?
list.pv <- do.call(rbind, lapply(1:50, function(tt){
  pv=NULL
  for(ff in c("hclust", "nmf", "svd", "subsampling")){
    Z0t <- initializeZt(Y, K=3L, force=TRUE, flavor = ff)
    res <- positiveFusedLasso(list(Y), Z0t, lambda) ## with the correct Z,Y= minus normal subclone,sample
    pv <- rbind(pv,c(it=tt, flavor= ff, modelFitStats(res) ))
  }
  return(pv)
}))
ggplot(as.data.frame(list.pv))+geom_boxplot(aes(x=flavor, y=loss %>% as.character %>% as.numeric() %>% round(2)))+ylab("Loss")
ggplot(as.data.frame(list.pv))+geom_boxplot(aes(x=flavor, y=PVE %>% as.character %>% as.numeric() %>% round(2)))+ylab("PVE")
### SVD does not provide good results


## so the pb is to reconstruct a normal y w/o an explicit normal subclone 

## - - - - - - - - - - - - - - - - - - - - 
## non-totally trivial W and remove normal component
## - - - - - - - - - - - - - - - - - - - - 
lambda <- 1e-4
W <- rbind(c(1, 1, 0)/2,
           c(1, 0, 1)/2,
            c(0,  1,1)/2)
stopifnot(all(rowSums(W)==1))
E <- matrix(rnorm(nrow(W)*ncol(Z[-1,]), sd=0), nrow=nrow(W), ncol=ncol(Z))
Y <- W %*% Z[-1,] + E

## Test with no normal component
Z0t <- initializeZt(Y, K=nrow(Y), forceNormal=FALSE)
res <- positiveFusedLasso(list(Y), Z0, lambda) 
round(res@W, 2) ## Identity
round(t(res@S$Z), 2)  ## Y
round(res@E$Y, 2)
round(res@mu, 2) 
modelFitStats(res)

## Test with  normal component
Z0Nt <- initializeZt(Y, K=nrow(Y), forceNormal=TRUE)
resForceNorm <- positiveFusedLasso(list(Y), Z0Nt, lambda) 
round(resForceNorm@W, 2) ## each profile is each feature
round(t(resForceNorm@S$Z), 2)  ## last row is completly different of normal, this logical because there is no normal component
round(resForceNorm@E$Y, 2)
round(resForceNorm@mu, 2) ## not 0 where the wild things are
modelFitStats(resForceNorm)

## Test with  normal component
Z0Nt <- initializeZt(Y, K=nrow(Z), forceNormal=TRUE)
resForceNorm <- positiveFusedLasso(list(Y), Z0Nt, lambda) ##  WtW is not invertible: no solution for this combinaison of lambda
round(resForceNorm@W, 2)## Component normal set to 0 (logical there is no normal component in the model)
round(t(resForceNorm@S$Z), 2)  ## Z is equal to Y
round(resForceNorm@E$Y, 2)
round(resForceNorm@mu, 2) ##  0 where the wild things are
modelFitStats(resForceNorm)
### When we force a compoenent to be normal, it is possible to solve the model with p>n

### Doesn't work very well, due to  number n = p so the trivial solution for W is 1 for each archetype which corresponds to each profile (non indentifiable problem)

###
list.pv <- do.call(rbind, lapply(1:10, function(tt){
  pv=NULL
  for(ff in c("hclust", "nmf", "svd", "subsampling")){
    Z0t <- initializeZt(Y, K=3L, force=TRUE, flavor = ff)
    res <- positiveFusedLasso(list(Y), Z0t, lambda) ## with the correct Z,Y= minus normal subclone,sample
    pv <- rbind(pv,c(it=tt, flavor= ff, modelFitStats(res) ))
  }
  return(pv)
}))
ggplot(as.data.frame(list.pv))+geom_boxplot(aes(x=flavor, y=loss %>% as.character %>% as.numeric() %>% round(2)))+ylab("Loss")
ggplot(as.data.frame(list.pv))+geom_boxplot(aes(x=flavor, y=PVE %>% as.character %>% as.numeric() %>% round(2)))+ylab("PVE")
## Here Hclust is better and more stable ...
