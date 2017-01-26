simulateSubclones <- function(len, nbClones, nBkp) {
    interval <- 1:(len - 1)
    u <- numeric(0)
    minLength <- 100
    while (length(u) < nBkp) {
        j <- sample(x = interval, size = 1, replace = FALSE)
        u <- c(u, j)
        b.inf <- max(1, j - minLength)
        b.sup <- min(len, j + minLength)
        v <- b.inf:b.sup
        interval <- setdiff(interval, v)
    }
    
    e <- c(sort(sample(size=4, x=1:(length(u)-1), replace=FALSE)), length(u))
    s <- c(1,e+1)[-(length(e)+1)]
    bkpsByClones <- mapply(function(ss,ee){
        sort(u[ss:ee])
    }, s, e)
    o <- order(u)
    
    ### Regions
    regNames <- c("(1,1)","(0,1)","(1,2)","(0,2)")
    pattern <- "\\(([0-9]),([0-9])\\)"
    regAnnot <- data.frame(region = regNames, 
                           freq = rep(1/4,4), stringsAsFactors = FALSE)
    regAnnot$C1 <- as.numeric(gsub(pattern, "\\1", regNames))
    regAnnot$C2 <- as.numeric(gsub(pattern, "\\2", regNames))
    
    candidateRegions <- function(regName) {
        if (is.null(regName)) 
            return(regAnnot[, "region"])
        reg <- subset(regAnnot, region == regName)
        d1 <- regAnnot[, "C1"] - reg[, "C1"]
        d2 <- regAnnot[, "C2"] - reg[, "C2"]
        ww <- which((d1 & !d2) | (!d1 & d2))
        regAnnot[ww, "region"]
    }
    
    tmp <- matrix("(1,1)", nrow=nbClones, ncol=len)
    
    for (rr in 1:nrow(tmp)) {
        start <- c(1, bkpsByClones[[rr]]+1)
        end <- c(bkpsByClones[[rr]], len)
        for(bb in 1:length(start)){
            if(bb==1){
                reg=NULL
            }else{
                reg <- tmp[rr,start[bb-1]]
            }      
            candReg <- candidateRegions(reg)
            reg <- sample(size=1, candReg)
            tmp[rr,start[bb]:end[bb]] <- reg
        }
    }
    
    regionsByClones <- sapply(1:length(bkpsByClones), function(bb){
        bkp <- c(bkpsByClones[[bb]],len)
        tmp[bb,bkp]
    })
    
    dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFraction=1)
    dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFraction=0)
    subClones <- c3co::buildSubclones(len, dataAnnotTP, dataAnnotN,
                                nbClones, bkpsByClones,regionsByClones)
    subClones
}


fllat <- function(dat, nb.arch.grid) {
    YTCNtoSeg <- t(sapply(dat, function(cc) cc$tcn))
    Y <- log2(YTCNtoSeg)-1
    result.pve <- FLLat::FLLat.PVE(na.omit(t(Y)), J.seq=nb.arch.grid) 
    res <- sapply(nb.arch.grid, function(pp) {
        rr <- FLLat::FLLat.BIC(na.omit(t(Y)), J=pp)
        W <- rr$opt.FLLat$Theta
        Z <- rr$opt.FLLat$Beta   
        resFLLAT <- list(res=list(Z=Z, W=t(W), Y.hat=(Y=t(W)%*%t(Z))), PVE=result.pve[which(pp==nb.arch.grid)])
    })
    res
}
