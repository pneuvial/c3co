## simulate subclones
subClones <- simulateSubclones(len, nbClones, nBkp)

## simulate copy number profiles from subclones and weight matrices
weightMats <- list()
dats <- list()
resC1C2 <- listenv::listenv()
resTCN <- listenv::listenv()
resFLLAT <- listenv::listenv()
for (ss in 1:nbSimu) {
    ## weight matrix
    M <- getWeightMatrix(70, 20, nbClones, n)
    weightMats[[ss]] <- M
    
    ## simulated profiles
    dat <- mixSubclones(subClones=subClones, M)
    dats[[ss]] <- dat
    ## c3co
    resC1C2[[ss]] %<-% c3co(dat, parameters.grid=parameters.grid, stat="C1C2")
    resTCN[[ss]] %<-% c3co(dat, parameters.grid=parameters.grid, stat="TCN")
    ## FLlat
    resFLLAT[[ss]] %<-% fllat(dat, nb.arch.grid=p.list)
}

saveRDS(subClones, file.path(pathSubClones, sprintf("subclones.rds")))
saveRDS(weightMats, file.path(pathWeights, sprintf("weightsMat.rds")))
saveRDS(dats, file.path(pathDat, sprintf("simu.rds")))
saveRDS(as.list(resC1C2), file.path(pathRes, sprintf("results_C3CO_C1C2.rds")))
saveRDS(as.list(resTCN), file.path(pathRes, sprintf("results_C3CO_TCN.rds")))
saveRDS(as.list(resFLLAT), file.path(pathRes, sprintf("results_FLLAT_TCN.rds")))

