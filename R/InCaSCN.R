#' Cancer subclone Inference function
#'
#' @param dat A list of data frame for each patient containing the total copy number \code{tcn}, the mirrored B allele fraction \code{dh}.
#' @param lambda1.grid A grid of real numbers which is the penalty coefficients for the fused lasso on the minor copy number dimension
#' @param lambda2.grid A grid of real numbers which is the penalty coefficients for the fused lasso on the major copy number dimension
#' @param nb.arch A vector of integers which is the number of archetypes in the model
#' @param stat TCN or C1C2
#' @param segment By defaut \code{TRUE} (segment data before inferring features)
#' @param output.dir directory to save segmentation and feature data
#' @param forceSeg \code{TRUE} of \code{FALSE} by default \code{FALSE} if \code{fileSeg="segData.rds"} already exists in \code{output.dir}
#' @param forceInferrence \code{TRUE} of \code{FALSE} by default \code{FALSE} if \code{fileFeat=featureData,p=p.rds} already exists in \code{output.dir}
#' @param init.random \code{TRUE} or \code{FALSE} by defaut \code{FALSE}. Initialization done by clustering
#' @return A list of archetypes (\code{Z} the total copy number matrix,\code{Z1} the minor copy number matrix and \code{Z2} the major copy number matrix), matrix weight \code{W} and the reconstructed minor and major copy numbers.
#' @examples
#' dataAnnotTP <- loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 3
#' bkps <- list(c(100,250)*10, c(150,400)*10,c(150,400)*10)
#' regions <-list(c("(0,3)", "(0,2)","(1,2)"), c("(1,1)", "(0,1)","(1,1)"), c("(0,2)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
#' M <- getWeightMatrix(100,0, 3, 15, sparse.coeff=0.7, contam.coeff=0.6, contam.max=2)
#' dat <- apply(M, 1, mixSubclones, subClones=datSubClone, fracN=NULL)
#' l1 <- seq(from=1e-6, to=1e-5, length=3)
#' l2 <- seq(from=1e-6, to=1e-5, length=3)
#' casResC1C2 <- InCaSCN(dat, lambda1.grid=l1, lambda2.grid=l2, nb.arch.grid=2:6)
#' casRes <- InCaSCN(dat, stat="TCN", lambda1.grid=l1, lambda2.grid=l2, nb.arch.grid=2:6)
#' @export
InCaSCN <- function(dat, lambda1.grid=NULL, lambda2.grid=NULL, nb.arch.grid=2:(length(dat)-1), stat="C1C2",output.dir="resultsInCaSCN", segment=TRUE,forceSeg=FALSE,forceInferrence=FALSE, init.random=FALSE){
  if(is.null(lambda1.grid)){
    lambda1.grid <- seq(from=1e-6, to=1e-5, length=10)
  }
  if(is.null(lambda2.grid)){
    lambda2.grid <- seq(from=1e-6, to=1e-5, length=10)
  }
  output.dir <- Arguments$getWritablePath(output.dir)
  if(segment){
    fileSeg <- file.path(output.dir,"segDat.rds")
    if(!file.exists(fileSeg)|| forceSeg){
      resSegmentation <- segmentData(dat, stat=stat)
      saveRDS(resSegmentation, file=fileSeg)
    }else{
      resSegmentation <- readRDS(file=fileSeg)
    }
    bkpList <- resSegmentation$bkp
  }else{
    cat("No segmentation: algorithm could take time\n")
  }
                                        # pca.res <- PCA(resSegmentation$Y, graph=FALSE)
                                        #  pp <- min(min(which(diff(pca.res$eig[[3]])<1e-3)), round(n/2))

  reslist <- list()
  BICp <- 1e8
  ##print(sprintf("nb.arch=%s",pp))
  cond <- TRUE
  it <- 1
  pp <- nb.arch.grid[it]
  while(cond){
    print(pp)
    BICp <- 1e8
    fileFeat <- file.path(output.dir, sprintf("featureData,p=%s.rds", pp))
    if(!file.exists(fileFeat)||forceInferrence){
      for(l1 in lambda1.grid){
        if(stat=="C1C2"){
          if(segment){
            Y1 <- t(resSegmentation$Y1)
            Y2 <- t(resSegmentation$Y2)
          }else{
            YTCNtoSeg <- t(sapply(dat, function(cc) cc$tcn))
            YDHtoSeg <- t(sapply(dat, function(cc) cc$dh))
            Y1 <- YTCNtoSeg*(1-YDHtoSeg)/2
            Y2 <- YTCNtoSeg*(1+YDHtoSeg)/2
          }
          for(l2 in lambda2.grid){
            n <- ncol(Y1)
            res <- InCaSCN:::positive.fused(Y1,Y2, pp, lambda1 = l1, lambda2 = l2,init.random)
            loss <- sum(((Y1+Y2)-(res$Y.hat$Y1+res$Y.hat$Y2))^2)
            kZ <- sum(apply(res$Z, 2, diff)!=0)
            BIC <-  n*ncol(Y1)*log(loss/(n*ncol(Y1)))+kZ*log(n*ncol(Y1))
            PVE <- 1-loss/(sum(((Y1+Y2)-rowMeans(Y1+Y2))^2))
            if(BIC<BICp){
              res.l <- list(BIC=BIC, PVE=PVE, res=res, param=list(nb.arch=pp, lambda1=l1, lambda2=l2), bkp=bkpList)
              BICp <- BIC
            }
          }
        }else{
          if(segment){
            Y <- t(resSegmentation$Y)

          }else{
            Y <- t(sapply(dat, function(cc) cc$tcn))
          }
          n <- ncol(Y)
          res <- InCaSCN:::positive.fused(Y, Y2=NULL, nb.arch=pp, lambda1 = l1,init.random)
          loss <- sum((Y-(res$Y.hat$Y1))^2)
          kZ <- sum(apply(res$Z, 2, diff)!=0)
          BIC <-  n*ncol(Y)*log(loss/(n*ncol(Y)))+kZ*log(n*ncol(Y))
          PVE <- 1-loss/(sum((Y-rowMeans(Y))^2))
          if(BIC<BICp){
            res.l <- list(BIC=BIC, PVE=PVE, res=res, param=list(nb.arch=pp, lambda1=l1), bkp=bkpList)
            BICp <- BIC
          }
        }
      }
      saveRDS(res.l, fileFeat)
    }else{
      res.l <- readRDS(fileFeat)
    }
    reslist[[it]] <- res.l
    it <- it+1
    pp <- nb.arch.grid[it]
    cond <- (sum(apply(res.l$res$Z, 2,function(ww) (sum(ww^2)<1e-3)))==0&sum(apply(res.l$res$W, 2,function(ww) (sum(ww^2)<1e-3)))==0&pp<max(nb.arch.grid))
  }
  return(reslist)
}
