context("Construction of weight matrices")

for (ii in 1:20) {
  M <- getWeightMatrix(prop.max = 100, prob.min = 0,
                       nb.arch = 5L, nb.samp = 20L,
                       sparse.coeff = 0.7,
                       contam.coeff = 0.6, contam.max = 2)
  expect_true(all(rowSums(M) <= 100))
  
  M <- getWeightMatrix(prop.max = 70, prob.min = 20,
                       nb.arch = 10L, nb.samp = 20L,
                       sparse.coeff = 0.7,
                       contam.coeff = 0.9,
                       contam.max = 20)
  
  expect_true(all(rowSums(M) <= 100))
}
