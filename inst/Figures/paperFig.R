library("c3co")
library("ggplot2")

## - - - - - - - - - - - - - - - - - - - - - - - - 
## parameters of subclones and samples
## - - - - - - - - - - - - - - - - - - - - - - - - 
nbClones <- 2
nBkp <- 7
### Breakpoints in subclones
bkpsByClones <-list(
    c(800, 1200, 2000, 2200), 
    c(400, 1800, 2200))
stopifnot(length(bkpsByClones) == nbClones)  ## sanity check

### Regions
regionsByClones <- list(
    c("(1,1)", "(1,2)", "(1,1)", "(0,1)", "(1,1)"),
    c("(0,1)", "(1,1)", "(1,2)", "(1,1)"))
stopifnot(length(regionsByClones) == nbClones)  ## sanity check

dataAnnotTP <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFraction=1)
dataAnnotN <- acnr::loadCnRegionData(dataSet="GSE13372", tumorFraction=0)

len <- 800*3 ## 3 because roughly 1/3 of heterozygous SNPs
stopifnot(all(sapply(bkpsByClones, max) < len))  ## sanity check
subClones <- buildSubclones(len, dataAnnotTP, dataAnnotN, nbClones, bkpsByClones, regionsByClones)

### Simulate Samples
n <- 2
M <- matrix(c(60, 20, 60, 0), ncol=n, byrow=TRUE)
dat <- mixSubclones(subClones=subClones, M)

df.tcn <- data.frame(tcn=c(dat[[1]]$tcn, dat[[2]]$tcn), pos=dat[[1]]$pos, sample=factor(rep(1:2, each=len)))

trueBkps <- sort(unlist(bkpsByClones))
trueBkps1 <- sort(bkpsByClones[[1]])
trueSeg1 <- matrixStats::binMeans(y=dat[[1]]$tcn, x=1:len, bx=c(1, unique(trueBkps), len))
trueSeg2 <- matrixStats:: binMeans(y=dat[[2]]$tcn, x=1:len, bx=c(1, trueBkps1, len))
df.trueSeg <- data.frame(y=c(trueSeg1, trueSeg1[length(trueSeg1)],
                             trueSeg2,trueSeg2[length(trueSeg2)]),
                         x=c(1, unique(trueBkps), len, 1, trueBkps1, len), 
                         sample=factor(rep(1:2, times=c(length(trueSeg1)+1,length(trueSeg2)+1))))
bkpsEst <- jointseg::jointSeg(Y=cbind(dat[[1]]$tcn,dat[[2]]$tcn), K=15)
EstBkps <- bkpsEst$bestBkp

EstSeg1 <- matrixStats::binMeans(y=dat[[1]]$tcn, x=1:len, bx=c(1,EstBkps,len))
EstSeg2 <- matrixStats::binMeans(y=dat[[2]]$tcn, x=1:len, bx=c(1,EstBkps,len))
df.EstSeg <- data.frame(y=c(EstSeg1, EstSeg1[length(EstSeg1)], EstSeg2, EstSeg2[length(EstSeg2)]),
                        x=rep(c(1, EstBkps, len), 2), 
                        sample=factor(rep(1:2, times=c(length(EstSeg1)+1, length(EstSeg2)+1))))

gp <- ggplot(df.tcn) + 
    geom_point(aes(y=tcn, x=pos), colour="#000000", cex=0.2) +
    facet_grid(sample~.) +
    theme_bw() + 
    ylim(c(0.5,3.5)) +
    scale_colour_manual(values=c("#1D702DAA", "#FFB300AA")) + 
    geom_step(data=df.EstSeg, aes(y=y, x=x, colour=sample), lty=1, lwd=1.5) + 
    ylab("Total copy number") +
    xlab("Genome position") + 
    geom_vline(xintercept=EstBkps, lty=2, colour="#FF0000", lwd=0.7)
gp

figPath <- R.utils::Arguments$getWritablePath("fig")
filename <- "exampleSimulatedProfile.png"
pathname <- file.path(figPath, filename)
ggsave(gp, filename=pathname, width=5, height=2.5)
