##########################################################################
# Download GEO data set GSE47077
##########################################################################
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

