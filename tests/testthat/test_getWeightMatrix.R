library("c3co")

context("Construction of weight matrices")

for (ii in 1:20) {
    M <- getWeightMatrix(100, 0, 5, 20, sparse.coeff=0.7, contam.coeff=0.6, contam.max=2);
    expect_true(all(rowSums(M)<=100))
    M <- getWeightMatrix(70, 20, 10, 20, sparse.coeff=0.7, contam.coeff=0.9, contam.max=20);
    expect_true(all(rowSums(M)<=100))
}