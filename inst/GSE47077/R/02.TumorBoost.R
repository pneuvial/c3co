##########################################################################
# TumorBoost Part
##########################################################################

## See aroma project to get more information about this part

library(aroma.cn)
rootPath <- "callData"
rootPath <- Arguments$getReadablePath(rootPath)

genotypeTag <- "NGC"
dataSet <- c("GSE47077,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY")
gsN <- AromaUnitGenotypeCallSet$byName(dataSet, tags=genotypeTag, chipType="*")

rootPath <- "totalAndFracBData"
rootPath <- Arguments$getReadablePath(rootPath)
ds <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath)
idxRK29 <- grep("*.RK29*", getNames(ds))
ds <- ds[idxRK29]
# Extract the arrays with this name, which should be the normal
idxN <- match("GSM1144476_RK29-N", getNames(ds))
dsN <- ds[idxN]
# and the tumor
idxT <- grep("*.RK29-R", getNames(ds))
dsT <- ds[idxT]


dsList <- list(normal=dsN, tumor=dsT, callsN=gsN)

sampleNames <- getNames(dsList$normal)


dummy <- lapply(dsList, FUN=function(ds) print(ds[[1]]))

tbn <- sapply(1:length(dsT), function (dt){
  print(dsList$tumor[dt])
  t <- TumorBoostNormalization(dsT=dsList$tumor[dt], dsN=dsList$normal, gcN=dsList$callsN, tags=c("*", "NGC"))
  dsTN <- process(t, verbose=log)
})

ugp <- getAromaUgpFile(dsList$tumor)
# Identify SNPs only
platform <- getPlatform(ugp)
if (platform == "Affymetrix") {
  require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix")
  snpPattern <- "^SNP"
} else if (platform == "Illumina") {
  snpPattern <- "^rs[0-9]"
} else {
  throw("Unknown platform: ", platform)
}
unf <- getUnitNamesFile(ugp)
## Extract Allele B fractions
dataSet <- c("GSE47077,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY")
dsTCN <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath)
idxRK29 <- grep("*.RK29*", getNames(dsTCN))
dsTCN <- dsTCN[idxRK29]
dataTCN <- extractPSCNArray(dsTCN)

dataSet <- c("GSE47077,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY,TBN,NGC")
dsBAF<- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath)
idxRK29 <- grep("*.RK29*", getNames(dsBAF))
dsBAF <- dsBAF[idxRK29]
dataBAF <- extractMatrix(dsBAF)

## Genotype
geno <- extractMatrix(gsN)

##########################################################################
# Build data for c3co
##########################################################################
dat <- lapply(2:14, function(ii) {
  print(dimnames(dataTCN)[[3]][ii])
  tcn <- 2*dataTCN[,"total",ii,drop=TRUE]/dataTCN[,"total",1,drop=TRUE]
  betaTN <- dataBAF[,(ii-1),drop=TRUE]
  dh <- 2*abs(betaTN-1/2)
  isHet <- (geno==1)
  dh[!isHet] <- NA
  pos <- ugp[,2]
  chromosome <- ugp[,1]
  isNA <- is.na(pos)
  df <- data.frame(tcn=tcn[!isNA], dh=dh[!isNA], pos=pos[!isNA], chr=chromosome[!isNA])
  o <- order(df$chr, df$pos)
  return(df[o,])
})

##########################################################################
# Save data
##########################################################################
names(dat) <- dimnames(data)[[3]][2:14]
path <- Arguments$getWritablePath("data")
patient <- "RK29"
saveRDS(dat, sprintf("%s/dat-%s.rds", path ,patient))
