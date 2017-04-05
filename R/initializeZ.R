#' Initialization of the c3co model parameters
#'
#' @export
#' @param Y1 A matrix containing the segmented minor copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param Y2 A matrix containing the segmented major copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param nb.arch An integer which is the number of archetypes in the model
#' @param init.random if you want to use random initialization set parameter to TRUE
#' @param flavor Statistic used to perform 'hclust' initialization. Should be either "C1+C2", "C1", or "C2"
#' @param verbose A logical value indicating whether to print extra information. Defaults to FALSE
#' @return A list with two components: \describe{\item{Z1}{A  \code{L} x \code{p} matrix, the initial value for the \code{L} minor copy numbers of the \code{p} latent features}\item{Z2}{A  \code{L} x \code{p} matrix, the initial value for the \code{L} major copy numbers of the \code{p} latent features}}
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100,250)*10, c(150,400)*10,c(150,400)*10)
#' regions <-list(c("(0,3)", "(0,1)","(1,2)"), 
#' c("(1,1)", "(0,1)","(1,1)"), c("(0,2)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
#' M <- rbind(c(40, 30, 0), 
#'   c(0, 70, 15),
#'   c(10, 0, 35),
#'   c(15, 0, 0),
#'   c(0, 0, 0))
#' simu <- mixSubclones(subClones=datSubClone, M)
#' seg <- segmentData(simu)
#' res <- initializeZ(seg$Y1, seg$Y2, nb.arch=4)
#' resC <- initializeZ(seg$Y1+seg$Y2, nb.arch=4)
#'
#' @importFrom stats dist hclust cutree
initializeZ <- function(Y1, Y2=NULL, nb.arch=ncol(Y1), init.random=FALSE, flavor=c("C1+C2", "C1", "C2"), verbose=FALSE) {
    n <- nrow(Y1) # number of samples
    L <- ncol(Y1) # number of loci/segments
    stopifnot(nb.arch<=n)
    flavor <- match.arg(flavor)
    
    if (is.null(Y2)){
        Y <- Y1
    } else {
        stopifnot(nrow(Y2)==n)  ## sanity check
        stopifnot(ncol(Y2)==L)  ## sanity check
        Y <- switch(flavor, 
                    "C1+C2"= Y1 + Y2,
                    "C1"= Y1,
                    "C2"= Y2)
    }
    if (!init.random){
        ## hierarchical agglomerative clustering on Y
        dd <- dist(Y)
#        dd <- as.dist(1 - 1/2*(cor(t(Y1)) + cor(t(Y2))))
        hc <- hclust(dd, method="ward.D")
        cluster <- cutree(hc, k=nb.arch)
        if (verbose) message("Clustering in ", nb.arch, " groups: ", cluster)
        
        ## subclones are initialized as intra-cluster averages
        Z.init <- sapply(split(as.data.frame(Y), f=cluster), FUN=colMeans)
        Z1.init <- sapply(split(as.data.frame(Y1), f=cluster), FUN=colMeans)
        Z2.init <- sapply(split(as.data.frame(Y2), f=cluster), FUN=colMeans)
        if (is.null(Y2)){
            Z1.init <- Z.init
            Z2.init <- NULL
        }
    } else {
        idxs <- sample(1:n, size=nb.arch, replace=FALSE)
        Z.init <- t(Y[idxs, , drop=FALSE])
        Z1.init <- t(Y1[idxs, , drop=FALSE])
        Z2.init <- t(Y2[idxs, , drop=FALSE])
        if (is.null(Y2)){
            Z2.init <- NULL
        }
    }
    Z <- list(Z=Z.init, Z1=Z1.init, Z2=Z2.init)
    return(Z)
}
