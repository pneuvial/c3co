##########################################################################
# Run c3co
##########################################################################
## We provide the results of the algorithm and the segmentation so if the user only wants to reproduce figure, run only this file

library("c3co")
library("RColorBrewer")
library("ggplot2")

patient <- "RK29"

### Plot PVE
pvePlot(resC3co@fit, ylim=c(0.70,1))
bestp <- 6
dataBest <- resC3co@fit[[bestp]]

### Plot W matrix
figPath <- R.utils::Arguments$getWritablePath("fig-GSE47077-RK29")
filename <- sprintf("heatmap,GSE47077,patient=%s.pdf", patient)
pathname <- file.path(figPath, filename)
pdf(pathname, width=13, height=8)
Wplot(resC3co, idxBest=bestp, rownamesW=sprintf("R%s",1:nrow(dataBest@W)))
dev.off()

### Plot Latent profiles
pathCHR <-  system.file("inst", "GSE47077", "data","minMaxposByCHR.rds", package = "c3co")
minMaxPos <- readRDS(pathCHR)

lengthCHR <- sapply(resC3co@bkp, length)
chrs <- sapply(1:22, function(cc) rep(cc, times=lengthCHR[cc]))

start <- c(1, cumsum(lengthCHR) + 1)

stats <- c("TCN", "Minor", "Major")
ylims <- list("TCN"=c(1, 3), 
              "Minor"=c(0, 2), 
              "Major"= c(1, 3))
for (ch in c(9, 14)) {
    for (what in stats) {
        df <- createZdf(resC3co, minMaxPos, chromosomes=ch, var=what, idxBest=bestp)
        df$position <- df$position/1e6 ## scale to Mb
        
        filename <- sprintf("archetypes,GSE47077,patient=%s,%s,chr=%s.pdf", patient, what, ch)
        pathname <- file.path(figPath, filename)
        gArch <- Zplot(df, ylab=what)
        ggsave(gArch, filename=pathname, width=7, height=3.5)
    }
}
