##########################################################################
# Allele-specific CRMAv2
##########################################################################
## Load GSE47077 data set

library("GEOquery")
library("R.utils")

dataSet <- "GSE47077";
chipType <- "GenomeWideSNP_6";
rpath <- "rawData"

## WARNING: The next command will attempt to download a 2 Gb file!
baseDir <- file.path("geoData", dataSet)
filePaths <- getGEOSuppFiles(dataSet, baseDir=baseDir)

fl <- file.path(baseDir, "filelist.txt")
library(dplyr)
dat <- read.table(fl, header=TRUE, as.is=TRUE)
tar <- sprintf("%s_RAW.tar", dataSet)
fnames <- dat[[tar]]

## untar "GSE47077_raw.tar" only for patient RK29
gzfiles <- grep(".*RK29-*", fnames, value=TRUE)
gzfiles

## gunzip the .CEL.gz files and save them to "rawData/GSE47077/GenomeWideSNP_6/"
path <- Arguments$getWritablePath(file.path(rpath, dataSet, chipType))
untar(file.path("GSE47077", "GSE47077_raw.tar"), files=gzfiles, exdir=path)

lapply(file.path(path, list.files(path)), gunzip)  ## is this necessary?

## Pre-processing using the Aroma project, see http://www.aroma-project.org and 
## http://www.aroma-project.org/blocks/doCRMAv2/

library("aroma.affymetrix")
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
print(cdf)

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, cdf=cdf)
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
resR <- doASCRMAv2(dataSet, verbose=verbose, cdf=cdf)


