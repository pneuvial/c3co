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
    dat <- apply(M, 1, mixSubclones, subClones=subClones, fracN=NULL)
    dats[[ss]] <- dat
    
    ## c3co
    resC1C2[[ss]] %<-% c3co(dat, nb.arch.grid=p.list, stat="C1C2")
    resTCN[[ss]] %<-% c3co(dat, nb.arch.grid=p.list, stat="TCN")
    
    ## FLlat
    resFLLAT[[ss]] %<-% fllat(dat, nb.arch.grid=p.list)
}
