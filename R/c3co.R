#' Cancer subclone Inference function
#'
#' @param dat A list of data frame for each patient. Data frame containing 
#' \describe{
#'   \item{tcn}{Total copy number}
#'   \item{dh}{Mirrored B allele fraction}
#'   \item{pos}{Position on the genome}
#'   \item{chr}{Chromosome}
#'   }
#' @param lambda1.grid,lambda2.grid A grid of real numbers which is the penalty coefficients for the fused lasso on the minor and major copy number dimension
#' @param nb.arch.grid A vector of integers which is the number of archetypes in the model
#' @param stat TCN or C1C2
#' @param saveResults By default \code{FALSE}, set parameter to TRUE if you want to save results in directory \code{output.dir}
#' @param output.dir Name of directory to save segmentation and feature data if \code{saveResults} is set to \code{TRUE}
#' @param pathSeg By default \code{NULL} if it is not \code{NULL} assuming that \code{fileSeg="segData.rds"} already exists in \code{pathSeg} and load segmentation.
#' @param init.random \code{TRUE} or \code{FALSE} by defaut \code{FALSE}. Initialization done by clustering
#' @param new.getZ \code{TRUE} if you want to parallelize inferrence of Minor and Major copy numbers (TRUE by default)
#' @param verbose A logical value indicating whether to print extra information. Defaults to FALSE
#' @return An object of class [\code{\linkS4class{c3coFit}}]
#' @examples
#' dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100,250)*10, c(150,400)*10,c(150,400)*10)
#' regions <-list(c("(0,3)", "(0,2)","(1,2)"), 
#' c("(1,1)", "(0,1)","(1,1)"), c("(0,2)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
#' M <- getWeightMatrix(100,0, 3, 15, sparse.coeff=0.7, contam.coeff=0.6, contam.max=2)
#' dat <- mixSubclones(subClones=datSubClone, M)
#' l1 <- seq(from=1e-6, to=1e-5, length=3)
#' l2 <- seq(from=1e-6, to=1e-5, length=3)
#' casResC1C2 <- c3co(dat, lambda1.grid=l1, lambda2.grid=l2, nb.arch.grid=2:6)
#' casRes <- c3co(dat, stat="TCN", lambda1.grid=l1, lambda2.grid=l2, nb.arch.grid=2:6)
#' @export
c3co <- function(dat, lambda1.grid=NULL, lambda2.grid=NULL, nb.arch.grid=2:(length(dat)-1), stat="C1C2", saveResults = FALSE, pathSeg=NULL, output.dir="results_c3co", init.random=FALSE, new.getZ=TRUE, verbose=FALSE){
  ## Sanity check
  if(!is.null(dat)){
     if(!is.list(dat)){
        stop("Argument 'dat' should be a list ")
      }
    checkCols <- lapply(dat, function (dd) {
      coln <- colnames(dd)
      ecn <- c("tcn", "dh", "pos", "chr") ## expected
      mm <- match(ecn, coln)
      if (any(is.na(mm))) {
        str <- sprintf("('%s')", paste(ecn, collapse="','"))
        stop("Argument 'dat' should contain columns named", str)
      }
    })
    expectedL <- nrow(dat[[1]]) 
    checklength <- lapply(dat, function (dd) {
      nrowDD <- nrow(dd)
      if (nrowDD != expectedL) {
        stop("Each data.frame in 'dat' should have the same size")
      }
    })
  }
  if(stat=="TCN"){
    new.getZ<-FALSE
  }
  if(is.null(lambda1.grid)){
    lambda1.grid <- seq(from=1e-6, to=1e-4, length=10)
  }
  if(is.null(lambda2.grid)){
    lambda2.grid <- seq(from=1e-6, to=1e-4, length=10)
  }
  if(saveResults){
    if(dir.exists(output.dir)){
      stop(sprintf("%s already exists. Remove it or change name directory to save the results"))
    }else{
      output.dir <- R.utils::Arguments$getWritablePath(output.dir)
    }
  }
  if(!is.null(pathSeg)){
    resSegmentation <- readRDS(file.path(pathSeg, "segDat.rds"))
  }else{
   resSegmentation <- segmentData(dat, stat=stat, verbose=verbose)
  }
  bkpList <- resSegmentation$bkp
  
  reslist <- methods::new("c3coFit")
  reslist@bkp <- bkpList
  reslist@segDat <- list(Y1=resSegmentation$Y1,Y2=resSegmentation$Y2,Y=resSegmentation$Y )
  
  cond <- TRUE
  it <- 1
  pp <- nb.arch.grid[it]
  while(cond){
    print(sprintf("number of Features = %s",pp))
    BICp <- 1e8
    fileFeat <- file.path(output.dir, sprintf("featureData,p=%s.rds", pp))
    if(file.exists(fileFeat)){
      message(sprintf("%s already exists: skipping", fileFeat))
      res.l <- readRDS(fileFeat)
    }else {
      for(l1 in lambda1.grid){
        if(stat=="C1C2"){
          Y1 <- t(resSegmentation$Y1)
          Y2 <- t(resSegmentation$Y2)
          for(l2 in lambda2.grid){
            res <- positive.fused(Y1,Y2, pp, lambda1 = l1, lambda2 = l2,init.random, new.getZ)
            if(res@BIC<BICp){
              res.l <- res
              BICp <- res@BIC
            }
          }
        }else if(stat=="TCN"){
          Y <- t(resSegmentation$Y)
          res <- positive.fused(Y, Y2 = NULL, nb.arch=pp, lambda1 = l1,init.random)
          if(res@BIC<BICp){
            res.l <- res
            BICp <- res@BIC
          }
        }else{
          stop("'stat' needs to be 'TCN' or 'C1C2'")
        }
      }
      if(saveResults){
        fileFeat <- file.path(output.dir, sprintf("featureData,p=%s.rds", pp))
        saveRDS(res.l, fileFeat)
        if (verbose) {
            message(sprintf("Results have been saved to file: %s", fileFeat))
        }
      }
    }
    reslist@fit[[it]] <- res.l
    it <- it+1
    pp <- nb.arch.grid[it]
    ## Z doesn't change anymore
    c1 <- sum(apply(res.l@S$Z, 2,function(ww) (sum(ww^2)<1e-3)))==0
    ## W doesn't change anymore
    c2 <- sum(apply(res.l@W, 2,function(ww) (sum(ww^2)<1e-3)))==0
    ## pp reach max of grid
    c3 <-  !is.na(pp)
    cond <- (c1& c2& c3)
  }
  return(reslist)
}



