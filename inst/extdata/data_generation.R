library(CVXR)
library(Matrix)

set.seed(784859)

## model settings
n <- 10L   ## Number of samples
J <- 11L   ## Number of segments
K <- 4L    ## Number of subclones
lambda <- 0.01 # penalty levalue for solving in Zt

## subclones
Z <- matrix(1, nrow = K, ncol = J)
Z[2, 2]    <- 2
Z[3, 5:6]  <- 2
Z[4, 9:10] <- 2

## weigths
W <- t(rmultinom(n, size=K, prob=c(1/8,1/8,1/4,1/2))/4)
E <- matrix(rnorm(n*J, sd = 0.1), nrow = n, ncol = J)
Y <- W %*% Z + E

## solving on W with CVXR
W_hat     <- CVXR::Variable(n, K)
my_objective <- CVXR::Minimize(sum((Y - W_hat %*% Z)^2))
my_constraints <- list(W_hat >= 0, W_hat %*% cbind(rep(1,K)) == cbind(rep(1,n)))
problem   <- CVXR::Problem(my_objective, my_constraints)
W_ref     <- solve(problem)$getValue(W_hat)

## solving on Z with CVXR
Z_hat <- CVXR::Variable(K, J)
D <- Matrix::bandSparse(
    J - 1, J, k = 1:0,
    list(rep(1, J - 1), rep(-1, J))
  )
loss    <- CVXR::Minimize(sum((Y - W %*% Z_hat)^2) + lambda * sum(abs(D %*% t(Z_hat))))
problem <- CVXR::Problem(loss)
Z_ref   <- solve(problem)$getValue(Z_hat)

data4ConsistencyTests <- list(Y = Y, Z = Z, W = W, E = E, Z_hat_ref = Z_ref, W_hat_ref = W_ref, lambda = lambda)

## saving results for testing consistency of get.W, get.Zt and eventually positiveFusedLasso
saveRDS(data4ConsistencyTests, "data4ConsistencyTests.rds")
