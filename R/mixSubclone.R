#' Create a mixture of subclones
#'
#' @param subClones The subclones used to create the mixture ideally from the function \code{buildSubClones}.
#' @param W A matrix of weights in percentage [0,100] (without fraction of normal cells)
#' @return The mixture
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 2
#' bkps <- list(c(100,250)*10, c(150,400)*10)
#' regions <-list(c("(0,3)", "(0,2)","(1,2)"), c("(1,1)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
#' mixture <- mixSubclones(datSubClone, c(20,30))
#' @export
mixSubclones <- function(subClones, W) {
    
    ## Sanity check
    idxHom <- which(subClones[[1]]$genotype!="0.5")
    sc <- sapply(seq_len(length(subClones)-1), FUN = function(i) {
        # Test if genotype of i is equal to genotype of j
        genoI <- subClones[[i]]$genotype
        for(j in (i+1):length(subClones)) {
            genoJ <- subClones[[j]]$genotype  
            try(if (sum(genoI==genoJ)!=length(genoI)) stop(sprintf("genotypes are not the same for sublones %s and %s", i, j )))
        }
    })
    
    ## compute parental copy numbers
    c1t <- sapply(subClones, FUN=function(ss) {
        dh <- 2*abs(ss$baft-1/2)
        c <- ss$ct*(1-dh)/2
    })
    
    c1n <-  sapply(subClones, FUN=function(ss) {
        dh <- 2*abs(ss$bafn-1/2)
        c <- ss$cn*(1-dh)/2
    })
    c2t <-  sapply(subClones, FUN=function(ss) {
        dh <- 2*abs(ss$baft-1/2)
        c <- ss$ct*(1+dh)/2
    })
    c2n <-  sapply(subClones, FUN=function(ss) {
        dh <- 2*abs(ss$bafn-1/2)
        c <- ss$cn*(1+dh)/2
    })
    
    if (is.vector(W)) {W <- matrix(W, nrow = 1)}
    df.res <- apply(W, MARGIN=1L, FUN=function(weights) {
        
        if (sum(weights)>100) { stop("Fraction Tumor upper than 100, please check weight vector") }
        fracN <- 100-sum(weights)
        weights <- weights/100
        fracN <- fracN/100
        
        c1 <- rowSums(cbind(sapply(seq(along=subClones), FUN=function(ii) {
            weights[ii]*c1t[,ii]
        }),fracN*rowMeans(c1n)))
        c2 <- rowSums(cbind(sapply(seq(along=subClones), FUN=function(ii) {
            weights[ii]*c2t[,ii]
        }),fracN*rowMeans(c2n)))
        tcn <- c1+c2
        
        c1[idxHom] <- NA_real_
        c1[idxHom] <- NA_real_
        c2[idxHom] <- NA_real_
        c2[idxHom] <- NA_real_
        dh <- (c2-c1)/(c2+c1)
        dh[idxHom] <- NA_real_
        data.frame(c1=c1,c2=c2, tcn=tcn, dh=dh,genotype=subClones[[1]]$genotype, chr=rep(1, times=length(c1)), pos=seq_along(c1))
    })
    
    return(df.res)
}
