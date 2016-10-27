library("aroma.affymetrix")
library(aroma.cn)
library(aroma.affymetrix)
library(cluster)
library(matrixStats)
library(jointseg)
library(ggplot2)
library(UPDohE)
library("NbClust")
library(GEOquery)
library(InCaSCN)

patient <- "OV03-17"

### Check if normalized data already exists
output.dir <- Arguments$getWritablePath(sprintf("resultsInCaSCN-%s-v2",patient))
fileSeg <- file.path(output.dir,"segDat.rds")

### Load data after aroma process
rootPath <- "totalAndFracBData"
rootPath <- Arguments$getReadablePath(rootPath)
dataSet <- c("GSE40546,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY")
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="GenomeWideSNP_6", paths=rootPath, verbose=TRUE)

## Get Serie Matrix file 
gse <- getGEO('GSE40546',GSEMatrix=TRUE)
identity <- pData(phenoData(gse[[1]]))

## Get the study id for each columns
id <- gsub("study id: ","",identity$characteristics_ch1)
## Get id for each study
studies <- unique(id)
## Group by patient
idxBySamples <- sapply(studies, function(oo) which(id==oo))

## Identity of focused patient
idxPatient <- idxBySamples[[patient]]
## GSM numbers of patient
gsm <- identity$geo_accession[idxPatient]
des <- identity$description[idxPatient]
fullName <- paste(gsm, des, sep="_")
## Extract data of patient
data <- extractPSCNArray(ds[fullName])

### Normalize TCN data
pathData <- "dataToDeconv"
filename <- "meanBySamples.rds"
filepath <- file.path(pathData, filename)
if(file.exists(filepath)){
  meansBySamples <- readRDS(filepath)
}else{
  meansBySamples <- sapply(idxBySamples, function(oo) rowMeans(data[,"total",oo, drop=FALSE]))
  saveRDS(meansBySamples, file.path(filepath))
}
mm <- rowMeans(meansBySamples[,-which(colnames(meansBySamples)==patient)])
CNData <- 2*data[,"total",]/mm
## Rename colnames
colnames(CNData) <- gsm

### SNPs
units <- seq(length=nbrOfUnits(ds[[1]]))
ugp <- getAromaUgpFile(ds)
platform <- getPlatform(ugp)
unf <- getUnitNamesFile(ugp)
units23 <- getUnitsOnChromosome(ugp, 23)
is23 <- is.element(units, units23)
units24 <- getUnitsOnChromosome(ugp, 24)
is24 <- is.element(units, units24)

## Formating datCN for InCaSCN after normalization
datCN <- apply(CNData, 2, function(dd){
  print(str(dd))
  df <- do.call(rbind, lapply(1:22, function(chromosome){
    chrTag <- sprintf("Chr%s", chromosome)
    units <- getUnitsOnChromosome(ugp, chromosome=chromosome)    
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
## Formating datBAF for InCaSCN
datB <- apply(data[,"fracB",], 2, function(dd){
  print(str(dd))
  df <- do.call(rbind, lapply(1:22, function(chromosome){
    chrTag <- sprintf("Chr%s", chromosome)
    units <- getUnitsOnChromosome(ugp, chromosome=chromosome)
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
if(!file.exists(fileSeg)){
  ## segmentation on all profiles
  resInCaSCN <- InCaSCN:::segmentData(datCN, stat="TCN")
  nbSamples <- length(gsm)

  datC1C2 <- lapply(1:nbSamples, function(jj){
    print(gsm[jj])
    dd <- datB[[jj]]
    ddC <- datCN[[jj]]
    dfFull <- data.frame(ddC, baf=dd$baf)
    ## Find normal regions
    normalRegByCHR <- sapply(1:22, function (cc){
      print(sprintf("chromosome :%s", cc))
      df <- subset(dfFull, chr==cc)
      posBkp <- resInCaSCN$bkp[[cc]]
      res <- findNormalReg(df, posBkp)
      print(str(res))
      if(is.na(res)){
        dis <- list(propAA=NA, propBB=NA, deltaAA=NA, deltaBB=NA)
      }else{
        ## Estimation of parameters of the normal regions
        posNorm <- res$pos
        dis <- distribEstimation(df, posNorm=posNorm)
      }
      return(dis)
    })## fin cc
    ## Summary on all normal regions
    propAA <- median(unlist(normalRegByCHR["propAA",]), na.rm=TRUE)
    propBB <- median(unlist(normalRegByCHR["propBB",]), na.rm=TRUE)
    deltaAA <- median(unlist(normalRegByCHR["deltaAA",]), na.rm=TRUE)
    deltaBB <- median(unlist(normalRegByCHR["deltaBB",]), na.rm=TRUE)
    propAB <- 1-(propAA+propBB)

    dis <- list(propAA=propAA, propBB=propBB, deltaAA=deltaAA,deltaBB=deltaBB)

    binC1C2 <- do.call(rbind,lapply(1:22, function (cc){
      print(sprintf("chromosome :%s", cc))
      df <- subset(dfFull, chr==cc)
      posBkp <- resInCaSCN$bkp[[cc]]
      ## estimation of DOH by segments
      res <- estimDoH(df, posBkp,dis, linTrick=TRUE,eps=0.3)
      ret <- data.frame(Y1=res$C1,Y2=res$C2, DH=res$dh, chrom=cc, Y=res$C2+res$C1)
      return(ret)
    }))
    return(binC1C2)
  })

  ## Save data 
  datC1C2toSave <- list()
  datC1C2toSave$Y1 <- do.call(cbind, lapply(datC1C2, function (dd) dd[,"Y1"]))
  datC1C2toSave$Y2 <- do.call(cbind, lapply(datC1C2, function (dd) dd[,"Y2"]))
  datC1C2toSave$Y <- do.call(cbind, lapply(datC1C2, function (dd) dd[,"Y"]))
  datC1C2toSave$bkp <- resInCaSCN$bkp

  saveRDS(datC1C2toSave, file=fileSeg)
}else{
  cat("File already exists")
}
