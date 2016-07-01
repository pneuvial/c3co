library("aroma.affymetrix")
library(aroma.cn)
library(aroma.affymetrix)
library(cluster)
library(matrixStats)
rootPath <- "totalAndFracBData"
rootPath <- Arguments$getReadablePath(rootPath)
dataSet <- c("GSE40546,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY")
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="GenomeWideSNP_6", paths=rootPath, verbose=TRUE)
sampleNames <- getNames(ds)
data <- extractPSCNArray(dsTCN=ds)

sampleNames <-  gsub("_....-...", "", getNames(ds)) # Work with the first sample only
identity <- read.table("dataToDeconv/GSE40546_series_matrix.csv", header=FALSE,nrows=14, sep=";", as.is=TRUE)
identity <- identity[,-1]
id <- gsub("study id: ","",identity[10,])
IDS <- unique(id)
idxBySamples <- sapply(IDS, function(oo) which(id==oo))
## meansBySamples <- sapply(idxBySamples, function(oo) rowMeans(data[,"total",oo, drop=FALSE]))
##saveRDS(meansBySamples, file.path(pathData,"meanBySamples.rds"))
stud <- "OV03-17"
deconvMarko <- function(stud){
  idxSample1 <- which(id==stud)
  gsm <- identity[2,idxSample1]
  ii <- sapply(gsm, grep, sampleNames)
  ## rm(data)
  pathData <- Arguments$getWritablePath("dataToDeconv")
  if(file.exists( file.path(pathData,sprintf("%s.rds", stud)))){
    dataS <- readRDS(file.path(pathData,sprintf("%s.rds", stud)))
    meansBySamples <- readRDS(file.path(pathData,"meanBySamples.rds"))
  }else{
    dataS <- data[,,ii]
    saveRDS(dataS, file.path(pathData,sprintf("%s.rds", stud)))
  }
#######################################################
##### CN after normalization and BAF
#######################################################
  mm <- rowMeans(meansBySamples[,-which(colnames(meansBySamples)==stud)])
  CNData <- 2*dataS[,"total",]/mm
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
  output.dir <- Arguments$getWritablePath(sprintf("resultsInCaSCN-%s",stud))
  resInCaSCN <- InCaSCN(dat,lambda1.grid=c(1e-6,2e-6,1e-5,2e-5),lambda2.grid=c(1e-6,2e-6,1e-5,2e-5), output.dir=output.dir)
  pathArch <- Arguments$getWritablePath("archData")
  saveRDS(resInCaSCN, file=file.path(pathArch, sprintf("archDataGSE40546,patient=%s.rds", stud)))
}

deconvMarko(stud)
