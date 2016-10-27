library("aroma.affymetrix")
library(aroma.cn)
library(aroma.affymetrix)
library(cluster)
library(matrixStats)

source("R/01.loadData.R")
output.dir <- Arguments$getWritablePath(sprintf("resultsInCaSCN-%s-v2",patient))
resInCaSCN <- InCaSCN(datCN,lambda1.grid=c(1e-6,2e-6,1e-5,2e-5),lambda2.grid=c(1e-6,2e-6,1e-5,2e-5), output.dir=output.dir,nb.arch.grid = 2:7)

