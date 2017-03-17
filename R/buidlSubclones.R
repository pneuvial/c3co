#' Build several subclones
#'
#' @param len The number of probes in subclones.
#' @param dataAnnotTP A data frame containing the data used for the resampling. Need the total copy number in tumor \code{c}, \code{b} the B allele fraction in tumor, the type of alteration in \code{region} and the genotype in reference \code{genotype}
#' @param dataAnnotN A data frame containing the data used for the resampling. Need the total copy number in normal \code{c}, \code{b} the B allele fraction in normal, the type of alteration in \code{region} and the genotype in reference \code{genotype}
#' @param nbClones The number of subclones
#' @param bkps The list of breakpoints by subclones
#' @param regions The list of altered regions by subclones
#' @return The list of subclone.
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 2
#' bkps <- list(c(100,250)*10, c(150,400)*10)
#' regions <-list(c("(0,3)", "(0,2)","(1,2)"), c("(1,1)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
#' 
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 2
#' bkps <- list(c(100,250)*10, c(150,400)*10)
#' regions <-list(c("(0,1)", "(0,2)","(1,2)"), c("(1,1)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN,  nbClones, bkps, regions)
#' @export
buildSubclones <- function(len, dataAnnotTP, dataAnnotN, nbClones, bkps=list(), regions=list()){
    if (nbClones != length(regions)){ 
        stop("Argument 'nbClones' should match 'length(regions)'")
    }
    if(nbClones!=length(bkps)){ 
        stop("Argument 'nbClones' should match 'length(bkp)'")
    }
    if (is.factor(dataAnnotTP$region)){
        dataAnnotTP$region <- as.character(dataAnnotTP$region)
    }
    genotype <- NULL; rm(genotype);
    if (is.factor(dataAnnotTP$genotype)){
        dataAnnotTP$genotype <- as.character(dataAnnotTP$genotype)
    }
    ## NB here we are forcing 1/3 hets and 2/3 homs. Why? /PN, 20170318
    idxAB <- sort(sample(x=1:len, size=len/3))
    idxAA <- sort(sample(x=(1:len)[-idxAB], size=len/3) )
    idxBB <- sort((1:len)[-c(idxAB,idxAA)])
    
    keepCols <- c("c", "b", "genotype", "region")
    dataAnnot <- data.frame(dataAnnotTP[, keepCols],
                            dataAnnotN[, c("c","b")])
    colnames(dataAnnot) <- c("ct", "baft", "genotype", "region", "cn", "bafn")
    
    subClone <- lapply(1:nbClones, function(ii){
        bkpsHet <- sapply(bkps[[ii]],function(bb) max(which(idxAB<=bb)))
        ssAB <- jointseg::getCopyNumberDataByResampling(length=length(idxAB), regData=subset(dataAnnot,genotype==0.5) , bkp=bkpsHet, regions=regions[[ii]])$profile
        bkpsAA <- sapply(bkps[[ii]],function(bb) max(which(idxAA<=bb)))
        ssAA <- jointseg::getCopyNumberDataByResampling(length=length(idxAA), regData=subset(dataAnnot,genotype==0) , bkp=bkpsAA, regions=regions[[ii]])$profile
        bkpsBB <- sapply(bkps[[ii]],function(bb) max(which(idxBB<=bb)))
        ssBB <- jointseg::getCopyNumberDataByResampling(length=length(idxBB), regData=subset(dataAnnot,genotype==1) , bkp=bkpsBB, regions=regions[[ii]])$profile
        pos <- c(idxAB, idxAA, idxBB)
        
        ss <- rbind(ssAB, ssAA,ssBB)
        ss$pos <- pos
        o <- order(pos)
        ss[o,]
    })
}

