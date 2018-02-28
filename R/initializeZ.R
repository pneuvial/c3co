#' Initialization of the latent features of the c3co model
#'
#' @param Y1 A matrix containing the segmented minor copy number (\code{n}
#'   patients in row and \code{L} segments in columns)
#'
#' @param Y2 A matrix containing the segmented major copy number (\code{n}
#'   patients in row and \code{L} segments in columns)
#'
#' @param p An integer value, the number of latent features in the model.
#'   Defaults to \code{min(dim(Y1))}
#'
#' @param flavor A character value specifying how initialization is perfomed.
#'   Defaults to \code{"hclust"}. See Details
#'
#' @param stat Statistic used to perform initialization. Should be either
#'   "C1+C2", "C1", or "C2"
#'   
#' @param forceNormal A logical value indicating whether a normal component is forced in initialization
#'   Defaults to FALSE
#'   
#' @param verbose A logical value indicating whether to print extra information.
#'   Defaults to FALSE
#'
#' @details The latent features (LF) are inferred as follows according to the
#'   value of argument 'flavor':
#'
#'   If \code{flavor=="hclust"} (the default), the LF are centers of clusters
#'   derived by hierarchical agglomerative clustering on the Euclidean distance
#'   between the input copy number profiles, and using Ward linkage
#'   (\code{\link[stats]{hclust}}).
#'
#'   If \code{flavor=="nmf"}, the LF are the _coefficients_ of the non-negative
#'   matrix factorization in \code{p} of the input copy number profiles.
#'
#'   If \code{flavor=="svd"}, the LF are the first \code{p} right singular
#'   vectors of the singular value decomposition of the input copy number
#'   profiles. The flavor is not recommended as it may produce matrices
#'   with non-positive entries
#'
#'   If  \code{flavor=="archetypes"}, the LF are defined using archetypal
#'   analysis.
#'
#'   If  \code{flavor=="subsampling"}, the LF are chosen at random among existing
#'   profiles.
#'
#' @references Gaujoux R and Seoighe C (2010). A flexible R package for
#'   nonnegative matrix factorization. BMC Bioinformatics, 11(1), pp. 367.
#'
#' @references Cutler A and Breiman L. (1994) Archetypal analysis.
#'   Technometrics, 36(4):338-3474.
#'
#' @return A list with two components: \describe{ \item{Z1}{A \code{L} x
#'   \code{p} matrix, the initial value for the \code{L} minor copy numbers of
#'   the \code{p} latent features} \item{Z2}{A \code{L} x \code{p} matrix, the
#'   initial value for the \code{L} major copy numbers of the \code{p} latent
#'   features} }
#'
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,1)", "(1,2)"),
#'                 c("(1,1)", "(0,1)", "(1,1)"),
#'                 c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len, nbClones, bkps, regions, dataAnnotTP, dataAnnotN)
#' M <- rSparseWeightMatrix(12, nbClones, 0.90)
#' simu <- mixSubclones(subClones=datSubClone, M)
#' seg <- segmentData(simu)
#' Y1 <- t(seg$Y1)
#' Y2 <- t(seg$Y2)
#'
#' resH <- initializeZ(Y1, Y2, p=nbClones)  ## corresponds to flavor "hclust")
#' resNMF <- initializeZ(Y1, Y2, p=nbClones, flavor="nmf")
#' \dontrun{
#' ## often fails because of singularities:
#' resArch <- initializeZ(Y1, Y2, p=nbClones, flavor="archetypes")
#' }
#' resSVD <- initializeZ(Y1, Y2, p=nbClones, flavor="svd")
#' resC <- initializeZ(Y1, Y2, p=nbClones, flavor="subsampling")
#'
#' resNMF1 <- initializeZ(Y1, p=nbClones, flavor="nmf")
#'
#' @importFrom stats dist hclust cutree
#' @export
initializeZ <- function(Y1, Y2=NULL, p=min(dim(Y1)),
                        flavor=c("hclust", "nmf", "archetypes", "svd", "subsampling"),
                        stat=c("C1+C2", "C1", "C2"), forceNormal=FALSE, verbose=FALSE) {
    n <- nrow(Y1) # number of samples
    L <- ncol(Y1) # number of loci/segments
   # stopifnot(p <= n)
    flavor <- match.arg(flavor)
    stat <- match.arg(stat)
    if(forceNormal) p <- p-1L
    if (is.null(Y2)) {
        Y <- Y1
    } else {
        ## sanity checks
        stopifnot(nrow(Y2) == n, ncol(Y2) == L)
        Y <- switch(stat,
                    "C1+C2" = Y1 + Y2,
                    "C1" = Y1,
                    "C2" = Y2)
    }
    if (flavor=="hclust") {
        dd <- dist(Y)
        hc <- hclust(dd, method="ward.D")
        initHclust <- function(Y, p) {
            cluster <- cutree(hc, k=p)
            t(sapply(split(as.data.frame(Y), f=cluster), FUN=colMeans))
        }
    } else if (flavor=="subsampling") {
        idxs <- sample(1:n, replace = FALSE)
        initSub <- function(Y, p) {
            Y[idxs[1:p], , drop=FALSE]
        }
    }

    initZ <- switch(flavor,
                    "nmf"=initNMF,
                    "svd"=initSVD,
                    "archetypes"=initArchetypes,
                    "hclust"=initHclust,
                    "subsampling"=initSub)

    if(forceNormal) {
      res <- list(Z1=cbind(t(initZ(Y1, p)), 1))
    }else{
      res <- list(Z1=t(initZ(Y1, p)))
    }
    if (!is.null(Y2)) {
      res$Z2 <- t(initZ(Y2, p))
      if(forceNormal) {
        res$Z2 <- cbind(res$Z2, 1)
      }
    }
    res
}

#' @importFrom NMF nmf coef
initNMF <- function(Y, p) {
    fit <- nmf(Y, rank = p)
    coef(fit)  ## NB: W is 'basis(fit)'
}

initSVD <- function(Y, p) {
    fit <- svd(Y, nu = 0, nv = p)
    t(fit$v)
}

#' @importFrom archetypes archetypes
initArchetypes <- function(Y, p) {
    fit <- archetypes(Y, k = p)
    fit$archetypes  ## Note: W is 'fit$alphas'
}
