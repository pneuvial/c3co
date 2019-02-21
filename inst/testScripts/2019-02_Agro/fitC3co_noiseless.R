set.seed(1288769)    
J <- 5 # segments

# subclones
z0 <- rep(1L, K)
z1 <- z0
z1[2] <- 0
z2 <- z0
z2[4] <- 0
Z <- rbind(z0, z1, z2)
nbBkps <- ncol(Z)-1

# weights
W <- diag(rep(1, nrow(Z)))  ## identity weights
Y <- W %*% Z                ## noiseless

# non-penalized
params <- list(lambda = 0, nb.arch = 1:3)
res <- fitC3co(Y, parameters.grid = params)
sapply(res$fit, modelFitStats) # model fit stats
Zhats <- lapply(res$fit, function(fit) t(fit@Zt$Z))

Zhats[[1]] - colMeans(Z)  ## close to 0
Zhats[[2]]
Zhats[[3]] - Z            ## close to 0

# penalized
params <- list(lambda = 1e-2, nb.arch = 1:K)
res <- fitC3co(Y, parameters.grid = params)
sapply(res$fit, modelFitStats) # model fit stats

Zhats <- lapply(res$fit, function(fit) t(fit@Zt$Z))

Zhats[[1]] - colMeans(Z)  ## less close to 0
Zhats[[2]]                
Zhats[[3]] - Z            ## less close to 0

# - - - - - - - - - - - - 
# more samples
# - - - - - - - - - - - - 
n <- 10
K <- nrow(Z)

W <- rSparseWeightMatrix(nb.samp = n, nb.arch = K, sparse.coeff=0.5)
Y <- W %*% Z                ## noiseless
Y <- as.matrix(Y)

# non-penalized
params <- list(lambda = 0, nb.arch = 1:8)
res <- fitC3co(Y, parameters.grid = params)
pvePlot2(res$config$best)
# ! negative PVE for a single sublclone. pb due to centering?