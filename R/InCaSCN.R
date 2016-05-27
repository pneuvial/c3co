#' Cancer subclone Inference function
#'
#' @param dat A list of data frame for each patient containing the total copy number \code{tcn}, the mirrored B allele fraction \code{dh}.
#' @param lambda1.grid A grid of real numbers which is the penalty coefficients for the fused lasso on the minor copy number dimension
#' @param lambda2.grid A grid of real numbers which is the penalty coefficients for the fused lasso on the major copy number dimension
#' @param nb.arch A vector of integers which is the number of archetypes in the model
#' @param stat
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
#' casResC1C2 <- InCaSCN(dat)
#' casRes <- InCaSCN(dat, stat="TCN")
#' @export
InCaSCN <- function(dat, lambda1.grid=NULL, lambda2.grid=NULL, nb.arch.grid=2:(length(dat)-1), stat="C1C2"){
  if(is.null(lambda1.grid)){
    lambda1.grid <- seq(from=1e-6, to=1e-5, length=10)
      }
  if(is.null(lambda2.grid)){
    lambda2.grid <- seq(from=1e-6, to=1e-5, length=10)
  }
  n <- length(dat)
  resSegmentation <- segmentData(dat, stat=stat)
 # pca.res <- PCA(resSegmentation$Y, graph=FALSE)
#  pp <- min(min(which(diff(pca.res$eig[[3]])<1e-3)), round(n/2))
  bkpList <- resSegmentation$bkp
  ### To do : parallel on l1
  reslist <- list()
  BICp <- 1e8
  ##print(sprintf("nb.arch=%s",pp))
  cond <- TRUE
  it <- 1
  pp <- nb.arch.grid[it]
  while(cond){
    print(pp)
    BICp <- 1e8
    for(l1 in lambda1.grid){
      if(stat=="C1C2"){
        Y1 <- t(resSegmentation$Y1)
        Y2 <- t(resSegmentation$Y2)
        for(l2 in lambda2.grid){
          res <- InCaSCN:::positive.fused(Y1,Y2, pp, lambda1 = l1, lambda2 = l2)
          loss <- sum(((Y1+Y2)-(res$Y.hat$Y1+res$Y.hat$Y2))^2)
          kZ <- sum(apply(res$Z, 2, diff)!=0)
          BIC <-  n*pp*log(loss/(n*pp))+kZ*log(n*pp)
          PVE <- 1-loss/(sum(((Y1+Y2)-rowMeans(Y1+Y2))^2))
          if(BIC<BICp){
            res.l <- list(BIC=BIC, PVE=PVE, res=res, param=list(nb.arch=pp, lambda1=l1, lamda2=l2), bkp=bkpList)
            BICp <- BIC
          }
        }
      }else{
        Y <- t(resSegmentation$Y)
        res <- InCaSCN:::positive.fused(Y, Y2=NULL, nb.arch=pp, lambda1 = l1)
        loss <- sum((Y-(res$Y.hat$Y1))^2)
        kZ <- sum(apply(res$Z, 2, diff)!=0)
        BIC <-  n*pp*log(loss/(n*pp))+kZ*log(n*pp)
        PVE <- 1-loss/(sum((Y-rowMeans(Y))^2))
        if(BIC<BICp){
          res.l <- list(BIC=BIC, PVE=PVE, res=res, param=list(nb.arch=pp, lambda1=l1), bkp=bkpList)
          BICp <- BIC
        }
      }
    }
  reslist[[it]] <- res.l
  it <- it+1
  pp <- nb.arch.grid[it]
  cond <- (sum(apply(res.l$res$Z, 2,function(ww) (sum(ww^2)<1e-3)))==0&sum(apply(res.l$res$W, 2,function(ww) (sum(ww^2)<1e-3)))==0&pp<max(nb.arch.grid))
  }
  return(reslist)
}
