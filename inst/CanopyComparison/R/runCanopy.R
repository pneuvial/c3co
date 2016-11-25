library(Canopy)
library(InCaSCN)
library(RColorBrewer)

message("load data set toy from Canopy")
projectname = 'toy'
data(toy2)

Y1 <- t(toy2$Wm)
Y2 <- t(toy2$WM)
lambda <- 1e-5

message("Run InCaSCN\n")
K = 2:8
rC1C2 <- lapply(K, function (kk) positive.fused(Y1,Y2, kk,lambda1 = lambda, lambda2 = lambda, init.random=FALSE))

n <- nrow(Y1)
loss <- sapply(rC1C2, function (rr) sum(((Y1+Y2)-(rr$Y.hat$Y1+rr$Y.hat$Y2))^2))
kZ <- sapply(rC1C2, function (rr) sum(apply(rr$Z, 2, diff)!=0))
PVE <- 1-loss/(sum(((Y1+Y2)-rowMeans(Y1+Y2))^2))
plot(K,PVE, type="l", ylim=c(0,1))
kbest <- 4

bestRes <- rC1C2[[which(K==kbest)]]
savePath <- Arguments$getWritablePath("data-Canopy")
saveRDS(bestRes, file=file.path(savePath, "resInCaSCN.rds"))
bestRes <- readRDS(file.path(savePath, "resInCaSCN.rds"))

W <- round(bestRes$W,2)[10:1,]
res.clust = hclust(dist(W),method="ward.D")
col = colorRampPalette(brewer.pal(9, 'GnBu'))(100)

figPath <- Arguments$getWritablePath("Figures-Canopy")
filename <- "heatmap,toy,incascn"
pdf(sprintf("%s/%s.pdf", figPath, filename), width=13, height=8)
heatmap.3(W, Rowv=FALSE,dendrogram="col", col=col,scale="none", cexCol=1.5, cexRow=1.5,margins = c(5,10), cellnote=W, notecol="black", key = FALSE)
dev.off()

rC1C2k4 <- lapply(1:50,function(ss) positive.fused(Y1,Y2, 4,lambda1 = lambda, lambda2 = lambda, init.random=TRUE))
loss <- sapply(rC1C2k4, function (rr) sum(((Y1+Y2)-(rr$Y.hat$Y1+rr$Y.hat$Y2))^2))

bestK4 <- rC1C2k4[[which.min(loss)]]

WK4 <- round(bestK4$W,2)[10:1,]
res.clust = hclust(dist(WK4),method="ward.D")
col = colorRampPalette(brewer.pal(9, 'GnBu'))(100)

figPath <- Arguments$getWritablePath("Figures-Canopy")
filename <- "heatmap,toy,incascn-v2"
pdf(sprintf("%s/%s.pdf", figPath, filename), width=13, height=8)
heatmap.3(WK4, Rowv=FALSE,dendrogram="col", col=col,scale="none", cexCol=1.5, cexRow=1.5,margins = c(5,10), cellnote=WK4, notecol="black", key = FALSE)
dev.off()



message("Canopy can take a little bit time")

R = toy2$R; X = toy2$X; WM = toy2$WM; Wm = toy2$Wm
epsilonM = toy2$epsilonM; epsilonm = toy2$epsilonm; Y = toy2$Y
K = 3:6; numchain = 5
sampchain = canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM,
 epsilonm = epsilonm, C = NULL, Y = Y, K = K,
 numchain = numchain, simrun = 50000, writeskip = 200,
 projectname = projectname, cell.line = FALSE,
 plot.likelihood = TRUE)


burnin = 100
thin = 10
# If pdf = TRUE, a pdf will be generated.
bic = canopy.BIC(sampchain = sampchain, projectname = projectname, K = K,
                 numchain = numchain, burnin = burnin, thin = thin, pdf = TRUE)
optK = K[which.max(bic)]

post = canopy.post(sampchain = sampchain, projectname = projectname, K = K,
                   numchain = numchain, burnin = burnin, thin = thin, 
                   optK = optK, post.config.cutoff = 0.05)
samptreethin = post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik = post[[2]]   # likelihoods of trees in samptree
config = post[[3]]
config.summary = post[[4]]
print(config.summary)

config.i = config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood.\n')
output.tree = canopy.output(post, config.i, C=NULL)

filename <- "heatmap,toy,canopy"
pdf.name = sprintf("%s/%s.pdf", figPath, filename)
canopy.plottree(output.tree, pdf.name =pdf.name, pdf = TRUE)
## To do compare InCaSCN to Canopy

tree <- output.tree
saveRDS(tree, file=file.path(savePath, "resCanopy.rds"))
tree <- readRDS(file.path(savePath, "resCanopy.rds"))

filename <- "heatmap,toy,canopy"
pdf.name = sprintf("%s/%s.pdf", figPath, filename)
pdf(sprintf("%s/%s.pdf", figPath, filename), width=13, height=8)
nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow = TRUE), widths = c(3, 3, 0), heights = c(0.5, 1, 0), respect =TRUE)
par(mar = c(1, 30, 1,10))
K = ncol(tree$CM)
plot(tree, label.offset = 0.1, type = "cladogram", direction = "d", 
     show.tip.label = FALSE)
nodelabels()
tiplabels()
snaedge = rep(NA, nrow(tree$sna))
for (k in 1:nrow(tree$sna)) {
  snaedge[k] = intersect(which(tree$edge[, 1] == tree$sna[k, 
                                            2]), which(tree$edge[, 2] == tree$sna[k, 3]))
}
cnaedge = rep(NA, nrow(tree$cna))
for (k in 1:nrow(tree$cna)) {
  cnaedge[k] = intersect(which(tree$edge[, 1] == tree$cna[k, 
                                            2]), which(tree$edge[, 2] == tree$cna[k, 3]))
}
edge.label = sort(unique(c(snaedge, cnaedge)))
tiplabels("Normal", 1, adj = c(0.2, 1.5), frame = "n", cex = 1.2, 
          col = 4)
tiplabels(paste("Clone", 1:(K - 2), sep = ""), 2:(K - 1), 
          adj = c(0.5, 1.5), frame = "n", cex = 1.2, col = 4)
tiplabels(paste("Clone", (K - 1), sep = ""), K, adj = c(0.8, 
                                                    1.5), frame = "n", cex = 1.2, col = 4)
par(mar = c(1, 30, 0.5, 9.5))
P = round(tree$P,2)
image(1:nrow(P), 1:ncol(P), axes = FALSE, ylab = "", xlab = "", 
      P, breaks = 0:100/100, col = col)
axis(4, at = 1:ncol(P), colnames(P), cex.axis = 1.2, las = 1, 
     tick = FALSE)
abline(h = seq(0.5, ncol(P) + 0.5, 1), v = seq(0.5, nrow(P) + 
                                                   0.5, 1), col = "grey")
for (i in 1:nrow(P)) {
  for (j in 1:ncol(P)) {
    txt.temp <- sprintf("%0.2f", P[i, j])
    if (P[i, j] <= 0.05 | P[i, j] >= 0.95) {
      text(i, j, txt.temp, cex = 1, col = "white")
    }
    else {
      text(i, j, txt.temp, cex = 1)
    }
  }
}
dev.off()
