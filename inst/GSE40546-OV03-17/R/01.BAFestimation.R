library("aroma.affymetrix")
library(aroma.cn)
library(aroma.affymetrix)
library(cluster)
library(matrixStats)
library(jointseg)
library(ggplot2)
library("NbClust")
### Test on Markowetz
rootPath <- "totalAndFracBData"
rootPath <- Arguments$getReadablePath(rootPath)
dataSet <- c("GSE40546,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY")
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="GenomeWideSNP_6", paths=rootPath, verbose=TRUE)
sampleNames <- getNames(ds)

sampleNames2 <-  gsub("_....-...", "", getNames(ds)) # Work with the first sample only
identity <- read.table("dataToDeconv/GSE40546_series_matrix.csv", header=FALSE,nrows=14, sep=";", as.is=TRUE)

identity <- identity[,-1]
id <- gsub("study id: ","",identity[10,])
IDS <- unique(id)
idxBySamples <- sapply(IDS, function(oo) which(id==oo))
## meansBySamples <- sapply(idxBySamples, function(oo) rowMeans(data[,"total",oo, drop=FALSE]))
##saveRDS(meansBySamples, file.path(pathData,"meanBySamples.rds"))
stud <- "OV03-17"
bafEst <- function(stud){
  idxSample1 <- which(id==stud)
  gsm <- identity[2,idxSample1]
  ii <- sapply(gsm, grep, sampleNames)
  
  data <- extractPSCNArray(dsTCN=ds[ii])
  pathData <- Arguments$getWritablePath("dataToDeconv")
  meansBySamples <- readRDS(file.path(pathData,"meanBySamples.rds"))

  mm <- rowMeans(meansBySamples[,-which(colnames(meansBySamples)==stud)])
  CNData <- 2*data[,"total",]/mm
  colnames(CNData) <- gsm
  units <- seq(length=nbrOfUnits(ds[[1]]))
  ugp <- getAromaUgpFile(ds)
  units23 <- getUnitsOnChromosome(ugp, 23)
  is23 <- is.element(units, units23)
  units24 <- getUnitsOnChromosome(ugp, 24)
  is24 <- is.element(units, units24)

#################################################
  ## deconv segmentation on all profiles
#################################################
  library(InCaSCN)
  dat <- apply(CNData, 2, function(dd){
    print(str(dd))
    df <- do.call(rbind, lapply(1:22, function(chromosome){
      chrTag <- sprintf("Chr%s", chromosome)
      units <- getUnitsOnChromosome(ugp, chromosome=chromosome)
      platform <- getPlatform(ugp)
      unf <- getUnitNamesFile(ugp)
      unitNames <- getUnitNames(unf, units=units)
      ## Identify SNP units
      snpPattern <- "^SNP"
      keep <- (regexpr(snpPattern, unitNames) != -1)
      units <- units[keep]
      pos <- getPositions(ugp, units=units)
      o <- order(pos)
      df <- data.frame(tcn=dd[units][o])
      df$pos <- pos[o]
      df$chr <- chromosome
      return(df)
    }))
  })

  datB <- apply(data[,"fracB",], 2, function(dd){
    print(str(dd))
    df <- do.call(rbind, lapply(1:22, function(chromosome){
      chrTag <- sprintf("Chr%s", chromosome)
      units <- getUnitsOnChromosome(ugp, chromosome=chromosome)
      platform <- getPlatform(ugp)
      unf <- getUnitNamesFile(ugp)
      unitNames <- getUnitNames(unf, units=units)
      ## Identify SNP units
      snpPattern <- "^SNP"
      keep <- (regexpr(snpPattern, unitNames) != -1)
      units <- units[keep]
      pos <- getPositions(ugp, units=units)
      o <- order(pos)
      df <- data.frame(baf=dd[units][o])
      df$pos <- pos[o]
      df$chr <- chromosome
      return(df)
    }))
  })

  resInCaSCN <- InCaSCN:::segmentData(dat, stat="TCN")
  distribNormBySample <- lapply(seq(along=datB), function(jj){
    print(jj)
    dd <- datB[[jj]]
    ddC <- dat[[jj]]
    normalRegByCHR <- sapply(1:22, function (cc){
      print(sprintf("chromosome :%s", cc))
      subDD <- subset(dd, chr==cc)
      subDDC <- subset(ddC, chr==cc)
      start <- c(1,ceiling(resInCaSCN$bkp[[cc]]))
      end <- c(floor(resInCaSCN$bkp[[cc]]), Inf)
      clust <- lapply(1:length(start), function (ii){
        ss <- which(subDD$pos>=start[ii]&subDD$pos<=end[ii])
        bsub <- as.matrix(na.omit(subDD[ss,"baf"]))
        csub <- subDDC[ss,"tcn"]
        if(length(bsub)>3){
          mc <-  tryCatch({
            NbClust(bsub, distance = "euclidean", min.nc = 2, max.nc = 5,method = "kmeans", index = "duda")
          },error=function(cond) {
            Best.nc <- c(Number_clusters=3, Value_Index=NA)
            Best.partition <- kmeans(bsub,  Best.nc["Number_clusters",])$cluster
            return(list(Best.nc, Best.partition))
          })
          nbClust <- mc$Best.nc
          part <- mc$Best.partition
        }else{
          print("clust=1")
          nbClust <- c(Number_clusters=1, Value_Index=NA)
          print(nbClust)
          part <- rep(1, length(bsub))
        }
        return(list(nbclust=nbClust,part=part, meanC=mean(csub)))
      })## fin seg
### Fin normal region
      meansC <- sapply(clust, function(dd) dd$meanC)
      part <- sapply(clust, function(dd) dd$part)
      prop <- sapply(clust, function(dd) table(dd$part)/sum(table(dd$part)))
      nbClust <- sapply(clust, function(dd) dd$nbclust["Number_clusters"])
      idxNorm <- which(nbClust==3)
      id3 <- sapply(prop[idxNorm], function(pp) sum((pp-1/3)^2))
      id2 <- meansC[idxNorm]
      if(length(rank(id3)+rank(id2))>0){
        norm <- idxNorm[which.min(rank(id3)+rank(id2))]
        ss <- which(subDD$pos>=start[norm]&subDD$pos<=end[norm])
        meanNormByGeno <- by(na.omit(subDD[ss,"baf"]), part[[norm]], mean)
### density of DH in normal region
        meanNormDH <- 2*abs(meanNormByGeno-1/2)
        normBB <- which(meanNormByGeno>0.80)
        normAA <- which(meanNormByGeno<0.20)
        propAA <- prop[norm][[1]][normAA]
        propBB <- prop[norm][[1]][normBB]
        res <- list(propAA=propAA, propBB=propBB, meanNormDHAA=meanNormDH[normAA],meanNormDHBB=meanNormDH[normBB],start=start[norm],end=end[norm])
        
      }else{
        cat("no normal region in this chromosome\n")
        res <- NA
      }
      return(res)
    })## fin cc
    return(normalRegByCHR)
  })



  propAA <- sapply(distribNormBySample, function (dis){
    sapply(dis, function (ch){
      if(!is.na(ch)){
        return(ch$propAA)
      }else{
        return(NA)
      }
    })
  })
  propBB <- sapply(distribNormBySample, function (dis){
    sapply(dis, function (ch){
      if(!is.na(ch)){
        return(ch$propBB)
      }else{
        return(NA)
      }
    })
  })

  meanNormDHAA <- sapply(distribNormBySample, function (dis){
    sapply(dis, function (ch){
      if(!is.na(ch)){
        return(ch$meanNormDHAA)
      }else{
        return(NA)
      }
    })
  })

  meanNormDHBB <- sapply(distribNormBySample, function (dis){
    sapply(dis, function (ch){
      if(!is.na(ch)){
        return(ch$meanNormDHBB)
      }else{
        return(NA)
      }
    })
  })
### Two regions are normal in all samples!

  start <- sapply(distribNormBySample, function (dis){
    sapply(dis, function (ch){
      if(!is.na(ch)){
        return(ch$start)
      }else{
        return(NA)
      }
    })
  })


  ComputeDHMeanHomInNorm <- mean(sapply(1:ncol(meanNormDHBB), function (ind){
    meanNormDHBB[,ind]*propBB[,ind]
  }), na.rm=TRUE)+ mean(sapply(1:ncol(meanNormDHAA), function (ind){
    meanNormDHAA[,ind]*propAA[,ind]
  }),na.rm=TRUE)

  propAB <- 1-(mean(sapply(1:ncol(meanNormDHBB), function (ind){
    propBB[,ind]
  }), na.rm=TRUE)+ mean(sapply(1:ncol(meanNormDHAA), function (ind){
    propAA[,ind]
  }),na.rm=TRUE))


  str(datB[[1]])
  datDH <- lapply(seq(along=datB), function (jj){
    print(jj)
    dd <- datB[[jj]]
    ddC <- dat[[jj]]
    binC1C2byCHR <- t(do.call(cbind, lapply(1:22, function (cc){
      print(sprintf("chromosome :%s", cc))
      subDD <- subset(dd, chr==cc)
      subDDC <- subset(ddC, chr==cc)
      start <- c(1,ceiling(resInCaSCN$bkp[[cc]]))
      end <- c(floor(resInCaSCN$bkp[[cc]]), Inf)
      binC1C2 <- sapply(1:length(start), function (ii){
        ss <- which(subDD$pos>=start[ii]&subDD$pos<=end[ii])
        dhsub <- (mean(2*abs(subDD[ss,"baf"]-1/2), na.rm=TRUE)-ComputeDHMeanHomInNorm)/propAB
        csub <- mean(subDDC[ss,"tcn"])
        c1 <- csub*(1-dhsub)/2
        c2 <- csub*(1+dhsub)/2
        if(c1<0) {c1=0}
        if(c2<c1){c2=c1} 
        return(c(Y1=c1, Y2=c2, Y=c1+c2))
      })
      return(binC1C2)
    })))
    res <- cbind(Y1=binC1C2byCHR[,"Y1"],Y2=binC1C2byCHR[,"Y"],Y=binC1C2byCHR[,"Y"])
    return(res)
  })
  datC1C2toSave <- list()
  datC1C2toSave$Y1 <- do.call(cbind, lapply(datDH, function (dd) dd[,"Y1"]))
  datC1C2toSave$Y2 <- do.call(cbind, lapply(datDH, function (dd) dd[,"Y2"]))
  datC1C2toSave$Y <- do.call(cbind, lapply(datDH, function (dd) dd[,"Y"]))
  datC1C2toSave$bkp <- resInCaSCN$bkp
  output.dir <- Arguments$getWritablePath(sprintf("resultsInCaSCN-%s",stud))
  fileSeg <- file.path(output.dir,"segDat.rds")
  saveRDS(datC1C2toSave, file=fileSeg)
}
bafEst(stud)
