context("Construction of weight matrices")
for (ii in 1:20) {
  W <- rSparseWeightMatrix(nb.samp = 20L, nb.arch = 5L, sparse.coeff = 0.7)
  expect_true(all(Matrix::rowSums(W, sparseResult = FALSE) <= 1))

  W <- rSparseWeightMatrix(nb.samp = 20L, nb.arch = 10L)
  expect_true(all(Matrix::rowSums(W, sparseResult = FALSE) <= 1))
}
