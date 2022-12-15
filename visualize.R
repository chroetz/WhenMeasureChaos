set.seed(1)

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


i0 <- 1#sample.int(n, 1)
x0 <- u[i0,]
#y0 <- u[i0,] + rnorm(3, sd=sd)
knn <- RANN::nn2(uSub[1:(nrow(uSub)-frames),], matrix(x0, nrow=1), k=k)

# par(mar=c(0,0,0,0), bg = 'lightgrey')
# plot(
#   NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
#   xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
#   ylab=NA, xlab=NA)
# text(x = mean(range(uSub[,dm[1]])), y = 40, "t = 0")
# points(uSub[, dm], pch=".")
diffs <- uSub[knn$nn.idx, dm]-rep(x0[dm], each=k)
angles <- atan2(diffs[,2], diffs[,1])
dists <- sqrt(rowSums(diffs^2))
cols <- hsv((angles+pi)/(2*pi), dists/max(dists), 1)
# points(uSub[knn$nn.idx, dm], pch=".", col=cols, cex=2)

for (i in c(1033,1043)) {
  cat(i,",",sep="")
  png(
    filename = file.path("imgs", sprintf("frame%04d.png", i)),
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
  points(uSub[knn$nn.idx+i, dm], pch=16, col=cols)
  dev.off()
}

