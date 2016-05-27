#' Create a mixture of subclones
#'
#' @param subClones The subclones used to create the mixture ideally from the function \code{buildSubClones}.
#' @param weights A vector of weights in percentage
#' @return The mixture
#' @examples
#' dataAnnotTP <- loadCnRegionData(dataSet="GSE11976", tumorFrac=1)
#' dataAnnotN <- loadCnRegionData(dataSet="GSE11976", tumorFrac=0)
#' len <- 500*10
#' nbClones <- 2
#' bkps <- list(c(100,250)*10, c(150,400)*10)
#' regions <-list(c("(0,3)", "(0,2)","(1,2)"), c("(1,1)", "(0,1)","(1,1)"))
#' datSubClone <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkps, regions)
#' mixture <- mixSubclones(datSubClone, c(20,30))
#' @export
mixSubclones <- function(subClones, weights, fracN=NULL){
  if(sum(weights)>100){stop("Fraction Tumor upper than 100, please check weight vector")}
  if(is.null(fracN)){
    fracN <- 100-sum(weights)
  }
  if(sum(weights,fracN)!=100){stop("Fraction Tumor+Normal is not equal to 100")}

  weights <- weights/100
  fracN <- fracN/100

  ## compute parental copy numbers
  c1t <- sapply(subClones, function(ss){
    dh <- 2*abs(ss$baft-1/2)
    c <- ss$ct*(1-dh)/2
  })

  c1n <-  sapply(subClones, function(ss){
    dh <- 2*abs(ss$bafn-1/2)
    c <- ss$cn*(1-dh)/2
  })
  c2t <-  sapply(subClones, function(ss){
    dh <- 2*abs(ss$baft-1/2)
    c <- ss$ct*(1+dh)/2
  })
  c2n <-  sapply(subClones, function(ss){
    dh <- 2*abs(ss$bafn-1/2)
    c <- ss$cn*(1+dh)/2
  })

  c1 <- rowSums(cbind(sapply(seq(along=subClones), function(ii){
    weights[ii]*c1t[,ii]
  }),fracN*rowMeans(c1n)))
  c2 <- rowSums(cbind(sapply(seq(along=subClones), function(ii){
    weights[ii]*c2t[,ii]
  }),fracN*rowMeans(c2n)))
  tcn <- c1+c2

  idxHom <- which(subClones[[1]]$genotype!="0.5")
  c1[idxHom] <- NA
  c1[idxHom] <- NA
  c2[idxHom] <- NA
  c2[idxHom] <- NA
  dh <- (c2-c1)/(c2+c1)
  dh[idxHom] <- NA

  return(data.frame(c1=c1,c2=c2, tcn=tcn, dh=dh,genotype=subClones[[1]]$genotype, chr=rep(1, length(c1)), pos=1:length(c1)))
}
