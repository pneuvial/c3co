##########################################################################
# Allele-specific CRMAv2
##########################################################################
## Load GSE47077 data set

library(GEOquery)
library(R.utils)
filePaths = getGEOSuppFiles("GSE47077")
listF <- list.files("GSE47077")
## untar "GSE47077_raw.tar" only for patient RK29
untar(file.path("GSE47077", "GSE47077_raw.tar"), files=grep(".*RK29-*", listF$Name, value=TRUE), exdir=pathToSaveCEL)

## gunzip the .CEL.gz files and save them to "rawData/GSE47077/GenomeWideSNP_6/"
pathToSaveCEL <- Arguments$getWritablePath("rawData/GSE47077/GenomeWideSNP_6/")
lapply(file.path(pathToSaveCEL, list.files(pathToSaveCEL)), gunzip)

## See aroma project to get more information about this part

library("aroma.affymetrix")
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE47077";
chipType <- "GenomeWideSNP_6";

cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
print(cdf)

## See aroma documentation for details
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, cdf=cdf)
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
resR <- doASCRMAv2(dataSet, verbose=verbose, cdf=cdf)


