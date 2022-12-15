set.seed(3)

attrFile <- "attractorLorenz63.RDS"
attractor <- readRDS(attrFile)
u <- attractor$u
n <- attractor$n

sd <- 1
k <- 20000
m <- 10
dm <- c(1,3)
subI <- seq(1, n, by = m)
uSub <- u[subI,]
frames <- 2000


i0 <- sample.int(n, 1)
x0 <- u[i0,]
knn <- RANN::nn2(uSub[1:(nrow(uSub)-frames),], matrix(x0, nrow=1), k=k)

idxs <- sort(knn$nn.idx)
dIdxs <- c(0, diff(idxs))
dIdxs[dIdxs != 1] <- 0
dIdxs <- 1-dIdxs
groups <- cumsum(dIdxs)
grTable <- table(groups)
gs <- (1:max(groups))[grTable >= 5]
groupMeanIdx <- sapply(gs, \(g) round(mean(idxs[groups == g])))
groupMeanVec <- uSub[groupMeanIdx, dm]
pca <- prcomp(groupMeanVec)
pcaDecomp <- (groupMeanVec-rep(pca$center, each=nrow(groupMeanVec))) %*% pca$rotation 
pcaValue <- pcaDecomp[,1]
pcaValue <- (pcaValue-min(pcaValue))/diff(range(pcaValue))*0.8

cols <- hsv(h = pcaValue, s=1, v=1)


for (i in 127:frames) {
#for (i in seq(0, frames, by=20)) {
  cat(i,",",sep="")
  png(
    filename = file.path("imgs2", sprintf("frame%04d.png", i)),
    width = 1024+512, height = 1024)
  par(mar = c(0,0,0,0), bg = 'lightgrey')
  plot(
    NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
    xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
    ylab=NA, xlab=NA)
  text(
    x = mean(range(uSub[,dm[1]])), y = 40,
    sprintf("t = %.2f", i/100),
    cex = 3)
  points(uSub[, dm], pch=".")
  for (l in seq_along(gs)) {
    lines(uSub[idxs[groups==gs[l]]+i, dm], col=cols[l], lwd=3)
  }
  dev.off()
}

