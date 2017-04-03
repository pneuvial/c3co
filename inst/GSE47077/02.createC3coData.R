##########################################################################
# Build input data for c3co
##########################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE47077"
patientID <- "RK29"

library("aroma.cn")  ## mostly for aroma.cn::callXXorXY
library("matrixStats")
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
rootPath <- "totalAndFracBData"
rootPath <- Arguments$getReadablePath(rootPath)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load AS-CRMAv2-normalized data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType="GenomeWideSNP_6", tags="ACC,ra,-XY,BPN,-XY,AVG,FLN,-XY")
patt <- sprintf(".*_%s-.*", patientID)
idx <- grep(patt, getNames(ds))
ds <- ds[idx]

sampleNames <- getNames(ds)
data <- extractPSCNArray(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call (naive) germline genotypes from normal BAF
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## chromosome/position annotation
ugp <- getAromaUgpFile(ds)
chr <- ugp[, 1, drop=TRUE]
pos <- ugp[, 2, drop=TRUE]

## Identify units on ChrX and ChrY
units23 <- getUnitsOnChromosome(ugp, 23)  ## identical to 'which(chr==23)'
units24 <- getUnitsOnChromosome(ugp, 24)  ## identical to 'which(chr==24)'

# Normal sample
idxN <- match("GSM1144476_RK29-N", sampleNames)
betaN <- data[,"fracB", idxN, drop=TRUE]

# Tumor sample
idxT <- grep("*.RK29-R", sampleNames)
betaT <- data[,"fracB", idxT, drop=TRUE]

## - - - - - - - - - - - - - -
## Germline genotypes
## - - - - - - - - - - - - - -
adjust <- 1

## Call gender
gender <- callXXorXY(betaN[units23], betaN[units24], adjust=adjust, from=0, to=1)
print(gender)  ## it's a male!
stopifnot(gender=="XY")  ## sanity check

## True total CN in normal sample (for a male)
cn <- rep(2, times=length(betaN))
cn[units23] <- 1
cn[units24] <- 1

mu <- callNaiveGenotypes(betaN, cn=cn, adjust=adjust, 
                          from=0, to=1, verbose=less(verbose,10))
print(table(mu, exclude=NULL))
print(table(mu, chr, exclude=NULL))

## - - - - - - - - - - - - - -
## Build input c3co data
## - - - - - - - - - - - - - -
totalN <- data[,"total", idxN,drop=TRUE]
totalT <- data[, "total", idxT]

keep <- which(!is.na(chr) & !is.na(pos))
chrK <- chr[keep]
posK <- pos[keep]
ord <- order(chrK, posK)

datList <- list()
for (ii in 1:ncol(betaT)) {
    ## Total copy number estimates
    tcn <- cn*totalT[, ii]/totalN  ## We are using 'cn' and not 2 so that X and Y CN estimates are centered around 1 (male), not 2!
    
    ## Allelic ratio estimates (after TumorBoost)
    betaTN <- aroma.light::normalizeTumorBoost(betaT[, ii], betaN, muN=mu)
    dh <- 2*abs(betaTN-1/2)
    isHet <- (mu==1/2)
    dh[!isHet] <- NA
    df <- data.frame(tcn=tcn[keep], dh=dh[keep], pos=posK, chr=chrK)
    datList[[ii]] <- df[ord, ]
}

## - - - - - - - - - - - - - -
## Save input data to file
## - - - - - - - - - - - - - -
names(datList) <- sampleNames[idxT]
path <- Arguments$getWritablePath("c3co-data")
filename <- sprintf("dat-%s.rds", patientID)
pathname <- file.path(path, filename)
saveRDS(datList, pathname)
