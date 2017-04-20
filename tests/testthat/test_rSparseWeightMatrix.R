context("Construction of weight matrices")
for (ii in 1:20) {
  M <- rSparseWeightMatrix(nb.samp = 20L,nb.arch = 5L, 
                        sparse.coeff = 0.7)
  expect_true(all(Matrix::rowSums(M,sparseResult = FALSE) <= 1))

  M <- rSparseWeightMatrix(nb.arch = 10, nb.samp = 20)
   
  expect_true(all(Matrix::rowSums(M,sparseResult = FALSE) <= 1))
}
