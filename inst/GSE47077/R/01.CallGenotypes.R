##########################################################################
# Call Genotype CRMAv2
##########################################################################

## See aroma project to get more information about this part
library("aroma.cn")
library(matrixStats)
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
rootPath <- "totalAndFracBData"
rootPath <- Arguments$getReadablePath(rootPath)
dataSet <- c("GSE47077,ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY")
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="GenomeWideSNP_6", paths=rootPath, verbose=TRUE)
idxRK29 <- grep("*.RK29*", getNames(ds))
ds <- ds[idxRK29]
sampleNames <- sort(unique(getNames(ds))) # Work with the first sample only
# Extract the arrays with this name, which should be the normal
idxN <- match("GSM1144476_RK29-N", getNames(ds))
dsN <- ds[idxN]
# and the tumor
idxT <- grep("*.RK29-R", getNames(ds))
dsT <- ds[idxT]


fullname <- paste(c(getFullName(dsN), "NGC"), collapse=",")
chipType <- getChipType(dsN, fullname=FALSE)
outPath <- file.path("callData", fullname, chipType)

units <- NULL
if (is.null(units)) {
  df <- dsN[[1]]
  units <- seq(length=nbrOfUnits(df))
}

# Identify units on ChrX and ChrY
ugp <- getAromaUgpFile(dsN)
units23 <- getUnitsOnChromosome(ugp, 23)
is23 <- is.element(units, units23)
units24 <- getUnitsOnChromosome(ugp, 24)
is24 <- is.element(units, units24)

kk <- 1
dfN <- dsN[[kk]]

tags <- getTags(dfN)
tags <- setdiff(tags, "fracB")
tags <- c(tags, "genotypes")
fullname <- paste(c(getName(dfN), tags), collapse=",")

filename <- sprintf("%s.acf", fullname)
gcPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE)

data <- extractPSCNArray(dsTCN=ds)
betaN <- data[units,"fracB","GSM1144476_RK29-N",drop=TRUE]

# Call genotypes
naValue <- as.double(NA)
fit <- NULL
mu <- rep(naValue, times=length(units))

# All but ChrX & ChrY in male
isDiploid <- (!(is23 | is24))
use <- which(isDiploid)
muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, from=0, to=1,
                                        verbose=less(verbose,10))
fit <- attr(muT, 'modelFit')
mu[use] <- muT
print(table(mu, exclude=NULL))

# Translate genotype calls in fracB space to (AA,AB,BB,...)
calls <- rep(as.character(NA), times=length(mu))
calls[mu ==   0] <- "AA"
calls[mu == 1/2] <- "AB"
calls[mu ==   1] <- "BB"
print(table(calls, exclude=NULL))

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Writing genotype calls (via temporary file)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- gcPathname
pathnameT <- sprintf("%s.tmp", pathname)
nbrOfUnits <- nbrOfUnits(dfN)
gfN <- AromaUnitGenotypeCallFile$allocate(pathnameT, platform=getPlatform(dfN), chipType=getChipType(dfN), nbrOfRows=nbrOfUnits)
footer <- readFooter(gfN)
footer$method <- "NaiveGenotypeCaller"
writeFooter(gfN, footer)

updateGenotypes(gfN, units=units, calls=calls)

res <- file.rename(pathnameT, pathname)
if (!isFile(pathname)) {
  throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname)
}
if (isFile(pathnameT)) {
  throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname)
}

gfN <- AromaUnitGenotypeCallFile(pathname)

print(gfN)
