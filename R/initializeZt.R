#' Initialization of the latent features of the c3co model
#'
#' @param Y1 A matrix containing the segmented minor copy number
#'   (n patients in rows and J segments in columns).
#'
#' @param Y2 A matrix containing the segmented major copy number
#'   (n patients in rows and J segments in columns).
#'
#' @param K An integer value, the number of latent features in the model.
#'   Defaults to `min(dim(Y1))`.
#'
#' @param flavor A character value specifying how initialization is performed.
#'   Defaults to [stats::hclust]. See Details.
#'
#' @param stat Statistic used to perform initialization. Should be either
#'   `"C1+C2"`, `"C1"`, or `"C2"`.
#'   
#' @param forceNormal A logical value indicating whether a normal component
#'   is forced in initialization.  Defaults to `FALSE`.
#'   
#' @param verbose A logical value indicating whether to print extra information.
#'   Defaults to `FALSE`.
#'
#' @details The latent features are inferred as follows according to the
#'   value of argument `flavor`:
#'
#'   If `flavor == "hclust"` (default), the latent features are centers of clusters
#'   derived by hierarchical agglomerative clustering on the Euclidean distance
#'   between the input copy number profiles, and using Ward linkage
#'   ([stats::hclust()]).
#'
#'   If `flavor == "nmf"`, the latent features are the _coefficients_ of the non-negative
#'   matrix factorization in `K` of the input copy number profiles.
#'
#'   If `flavor == "svd"`, the latent features are the first `K` right singular
#'   vectors of the singular value decomposition of the input copy number
#'   profiles. The flavor is not recommended as it may produce matrices
#'   with non-positive entries
#'
#'   If  `flavor == "archetypes"`, the latent features are defined using archetypal
#'   analysis.
#'
#'   If  `flavor == "subsampling"`, the latent features are chosen at random among existing
#'   profiles.
#'
#' @references Gaujoux R and Seoighe C (2010). A flexible R package for
#'   nonnegative matrix factorization. BMC Bioinformatics, 11(1), pp. 367.
#'
#' @references Cutler A and Breiman L. (1994) Archetypal analysis.
#'   Technometrics, 36(4):338-3474.
#'
#' @return A list with two components:
#' \describe{
#'   \item{`Z1`}{An J-by-`K` matrix, the initial value for the J minor
#'    copy numbers of the `K` latent features}
#'   \item{`Z2`}{An J-by-`K` matrix, the initial value for the J major
#'    copy numbers of the `K` latent features}
#' }
#' _Warning_: Note that the `Z1` and `Z2` components are actually the
#' transposed version of those matrices.  This notation mistake will be
#' fixed in a future release.
#'
#' @examples
#' ## Simulate locus-level (C1,C2) copy-number data
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10  ## Number of loci
#' K <- 3L        ## Number of subclones
#' n <- 12L       ## Number of samples
#' bkps <- list(c(100, 250)*10, c(150, 400)*10, c(150, 400)*10)
#' regions <- list(c("(0,3)", "(0,1)", "(1,2)"),
#'                 c("(1,1)", "(0,1)", "(1,1)"),
#'                 c("(0,2)", "(0,1)", "(1,1)"))
#' datSubClone <- buildSubclones(len=len, nbClones=K, bkps=bkps, regions=regions,
#'                               dataAnnotTP=dataAnnotTP, dataAnnotN=dataAnnotN)
#' W <- rSparseWeightMatrix(nb.samp=n, nb.arch=K, sparse.coeff=0.90)
#' simu <- mixSubclones(subClones=datSubClone, W=W)
#'
#' ## Segment the copy-number data
#' seg <- segmentData(simu)
#'
#' ## Initialize C3CO model
#' Y1 <- t(seg$Y1)
#' Y2 <- t(seg$Y2)
#'
#' resH <- initializeZt(Y1, Y2, K=K)  ## corresponds to flavor "hclust")
#' resNMF <- initializeZt(Y1, Y2, K=K, flavor="nmf")
#' \dontrun{
#' ## often fails because of singularities:
#' resArch <- initializeZt(Y1, Y2, K=K, flavor="archetypes")
#' }
#' resSVD <- initializeZt(Y1, Y2, K=K, flavor="svd")
#' resC <- initializeZt(Y1, Y2, K=K, flavor="subsampling")
#'
#' resNMF1 <- initializeZt(Y1, K=K, flavor="nmf")
#'
#' @importFrom stats dist hclust cutree
#' @export
initializeZt <- function(Y1, Y2=NULL, K=min(dim(Y1)),
                        flavor=c("hclust", "nmf", "archetypes", "svd", "subsampling"),
                        stat=c("C1+C2", "C1", "C2"), forceNormal=FALSE, verbose=FALSE) {
    n <- nrow(Y1) # number of samples
    J <- ncol(Y1) # number of loci/segments
    stop_if_not(is.numeric(K), length(K) == 1L, is.finite(K), K > 0L)
   # stop_if_not(K <= n)
    flavor <- match.arg(flavor)
    stat <- match.arg(stat)
    if (forceNormal) K <- K-1L
    if (is.null(Y2)) {
        Y <- Y1
    } else {
        ## sanity checks
        stop_if_not(nrow(Y2) == n, ncol(Y2) == J)
        Y <- switch(stat,
                    "C1+C2" = Y1 + Y2,
                    "C1" = Y1,
                    "C2" = Y2)
    }
    if (flavor == "hclust") {
        dd <- dist(Y)
        hc <- hclust(dd, method="ward.D")
        initHclust <- function(Y, K) {
            cluster <- cutree(hc, k=K)
            t(sapply(split(as.data.frame(Y), f=cluster), FUN=colMeans))
        }
    } else if (flavor == "subsampling") {
        idxs <- sample(1:n, replace = FALSE)
        initSub <- function(Y, K) {
            Y[idxs[1:K], , drop=FALSE]
        }
    }

    initZ <- switch(flavor,
                    "nmf"=initNMF,
                    "svd"=initSVD,
                    "archetypes"=initArchetypes,
                    "hclust"=initHclust,
                    "subsampling"=initSub)

    res <- list(Z1 = t(initZ(Y1, K = K)))
    if (!is.null(Y2)) res$Z2 <- t(initZ(Y2, K = K))

    if (forceNormal) {
      for (name in names(res)) res[[name]] <- cbind(res[[name]], 1)
    }
    
    ## Sanity checks
    for (name in names(res)) {
      stop_if_not(nrow(res[[name]]) == J, ncol(res[[name]]) == K)
    }
    
    res
}

#' @importFrom NMF nmf coef
initNMF <- function(Y, K) {
    fit <- nmf(Y, rank = K)
    coef(fit)  ## NB: W is 'basis(fit)'
}

initSVD <- function(Y, K) {
    fit <- svd(Y, nu = 0L, nv = K)
    t(fit$v)
}

#' @importFrom archetypes archetypes
initArchetypes <- function(Y, K) {
    fit <- archetypes(Y, k = K)
    fit$archetypes  ## Note: W is 'fit$alphas'
}
