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
# Extract the two arrays with this name, which should be the tumor and the normal
idxN <- match("GSM1144476_RK29-N", getNames(ds))
dsN <- ds[idxN]
dsT <- ds[-idxN]


dsList <- list(normal=dsN, tumor=dsT, callsN=gsN)

sampleNames <- getNames(dsList$normal)


dummy <- lapply(dsList, FUN=function(ds) print(ds[[1]]))

tbn <- sapply(1:length(dsT), function (dt){
  print(dsList$tumor[dt])
  t <- TumorBoostNormalization(dsT=dsList$tumor[dt], dsN=dsList$normal, gcN=dsList$callsN, tags=c("*", "NGC"))
  dsTN <- process(t, verbose=log)
}
              )

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
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="*", paths=rootPath)
data <- extractPSCNArray(ds)
#betaN <- data[units,"fracB","GSM1144476_RK29-N",drop=TRUE]
#betaT <- data[units,"fracB",-1,drop=TRUE]
#CN <- 2*data[units,"total",-1,drop=TRUE]/data[units,"total",1,drop=TRUE]


dataSet <- c("GSE47077,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY,TBN,NGC")
ds2 <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType="*", paths=rootPath)
data2 <- extractMatrix(ds2)
#betaTN <- data2[units,,drop=TRUE]

## Genotype
geno <- extractMatrix(gsN)

##########################################################################
# Build data for c3co
##########################################################################
dat <- lapply(2:14, function(ii) {
  print(dimnames(data)[[3]][ii])
  tcn <- 2*data[,"total",ii,drop=TRUE]/data[,"total",1,drop=TRUE]
  betaTN <- data2[,(ii-1),drop=TRUE]
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
