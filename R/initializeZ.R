#' Initialization of the c3co model parameters
#'
#' @export
#' @param Y1 A matrix containing the segmented minor copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param Y2 A matrix containing the segmented major copy number (\code{n} patients in row and \code{L} segments in columns)
#' @param nb.arch An integer which is the number of archetypes in the model
#' @param init.random if you want to use random initialization set parameter to TRUE
#' @param verbose A logical value indicating whether to print extra information. Defaults to FALSE
#' @return A list with two components: \describe{\item{Z1}{initial value for the minor copy number of the latent features}\item{Z2}{initial value for the major copy number of the latent features}}
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
#'
initializeZ <- function(Y1, Y2, nb.arch, init.random=FALSE, verbose=FALSE) {
    
    ## problem dimensions
    n <- nrow(Y1) # number of individuals
    L <- ncol(Y1) # number of loci/segments
    ## _______________________________________________________
    ## STEP 0: INITIALIZATION
    if (is.null(Y2)){
        Y <- Y1
    } else {
        Y <- Y1 + Y2
    }
    if (!init.random){
        ## initializing Z by clustering on Y
        hc <- stats::hclust(stats::dist(Y),method="ward.D")
        cluster <- stats::cutree(hc, nb.arch)
        ## averaging the Y over the clusters to initialize the archetypes
        Z1.init <- sapply(split(as.data.frame(Y1), cluster), colMeans)
        if (!is.null(Y2)){
            Z2.init <- sapply(split(as.data.frame(Y2), cluster), colMeans)
        } else {
            Z2.init <- NULL
        }
    } else {
        ii <- sample(1:n,nb.arch, replace=FALSE)
        Z1.init <- t(Y1[ii, ])
        Z2.init <- NULL
        if (!is.null(Y2)){
            Z2.init <- t(Y2[ii, ])
        }
    }
    Z <- list(Z1=Z1.init, Z2=Z2.init)
    return(Z)
}
