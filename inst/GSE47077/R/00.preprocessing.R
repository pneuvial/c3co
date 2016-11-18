##########################################################################
# Allele-specific CRMAv2
##########################################################################

## See aroma project to get more information about this part

library("aroma.affymetrix");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSet <- "GSE47077";
chipType <- "GenomeWideSNP_6";

cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full")
print(cdf)

## See aroma documentation for details
csR <- AffymetrixCelSet$byName(dataSet, chipType=chipType, cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AS-CRMAv2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
resR <- doASCRMAv2(dataSet, verbose=verbose, cdf=cdf);
dsC <- resR$total
print(dsC)


