##########################################################################
# Allele-specific CRMAv2
##########################################################################

## Pre-processing using the Aroma project, see http://www.aroma-project.org and 
## http://www.aroma-project.org/blocks/doCRMAv2/

dataSet <- "GSE47077";
chipType <- "GenomeWideSNP_6";


library("aroma.affymetrix")
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType(chipType, tags="Full")
print(cdf)

csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, cdf=cdf)
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
resR <- doASCRMAv2(dataSet, verbose=verbose, cdf=cdf)
