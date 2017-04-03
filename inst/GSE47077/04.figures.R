##########################################################################
# Reproducing the manuscript figures
##########################################################################

library("c3co")
library("RColorBrewer")
library("ggplot2")

dataSet <- "GSE47077"
patientID <- "RK29"

ptag <- sprintf("%s,patient=%s", dataSet, patientID)
figPath <- file.path("fig", ptag)
figPath <- R.utils::Arguments$getWritablePath(figPath)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot PVE
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
filename <- sprintf("PVE,%s.pdf", ptag)
pathname <- file.path(figPath, filename)
pdf(pathname)
pvePlot(resC3co@fit, ylim=c(0.70,1))
dev.off()


bestp <- 6
dataBest <- resC3co@fit[[bestp]]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot matrix W
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

filename <- sprintf("heatmap,%s.pdf", ptag)
pathname <- file.path(figPath, filename)
pdf(pathname, width=13, height=8)
Wplot(resC3co, idxBest=bestp, rownamesW=sprintf("R%s",1:nrow(dataBest@W)))
dev.off()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot latent profiles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathCHR <-  system.file("inst", dataSet, "data","minMaxposByCHR.rds", package = "c3co")
minMaxPos <- readRDS(pathCHR)

lengthCHR <- sapply(resC3co@bkp, length)
chrs <- sapply(1:22, function(cc) rep(cc, times=lengthCHR[cc]))
start <- c(1, cumsum(lengthCHR) + 1)
stats <- c("TCN", "Minor", "Major")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# selected chromosomes:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chs <- c(9, 14)
chrTag <- sprintf("chr=%s", paste(chs, collapse="-"))

df <- createZdf(resC3co, minMaxPos, chromosomes=chs, var=stats, idxBest=bestp)
df$position <- df$position/1e6 ## scale to Mb
filename <- sprintf("features,%s,%s.pdf", ptag, chrTag)
pathname <- file.path(figPath, filename)
gArch <- Zplot(df)
ggsave(gArch, filename=pathname, width=10, height=10)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# all chromosomes:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (chs in 1:24) {
    chrTag <- sprintf("chr=%s", paste(chs, collapse="-"))

    df <- createZdf(resC3co, minMaxPos, chromosomes=chs, var=stats, idxBest=bestp)
    df$position <- df$position/1e6 ## scale to Mb
    filename <- sprintf("features,%s,%s.pdf", ptag, chrTag)
    pathname <- file.path(figPath, filename)
    gArch <- Zplot(df)
    ggsave(gArch, filename=pathname, width=10, height=10)
}