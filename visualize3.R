
attrFile <- "attractorLorenz63.RDS"
attractor <- readRDS(attrFile)
u <- attractor$u
n <- attractor$n

sd <- 5
k <- 20000
m <- 10
dm <- c(1,3)
subI <- seq(1, n, by = m)
uSub <- u[subI,]
frames <- 2000

deltaI <- 100

cols <- c("#000000FF", viridis::inferno(1000))


png(
  filename = file.path("imgs3", "attractor.png"),
  width = 1024+512, height = 1024)
par(mar = c(0,0,0,0), bg = 'lightgrey')
plot(
  NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
  xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
  ylab=NA, xlab=NA)
points(uSub[, dm], pch=".")
dev.off()



for (seed in c(14:20)) {
    
  set.seed(seed)
  
  i0 <- sample.int(n, 1)
  x0 <- u[i0,]
  x1 <- u[i0+deltaI,]
  y0 <- x0 + rnorm(3, sd=sd)
  y1 <- x1 + rnorm(3, sd=sd)
  
  y0Dst <- DEEBesti:::distToVec(uSub, y0)
  y0i <- which.min(y0Dst)
  y0Proj <- uSub[y0i,]
  nSub <- nrow(uSub) - deltaI
  y0ProjDst <- DEEBesti:::distToVec(uSub[1:nSub,], y0Proj)
  y0W <- dnorm(y0ProjDst, sd=sd)
  y0ColI <- ceiling(y0W/max(y0W)*1000)+1
  
  y1Dst <- DEEBesti:::distToVec(uSub, y1)
  y1i <- which.min(y1Dst)
  y1Proj <- uSub[y1i,]
  y1ProjDst <- DEEBesti:::distToVec(uSub[(1:nSub)+deltaI,], y1Proj)
  y1W <- dnorm(y1ProjDst, sd=sd)
  y1ColI <- ceiling(y1W/max(y1W)*1000)+1
  
  y01W <- y0W*y1W
  y01ColI <- ceiling(y01W/max(y01W)*1000)+1
  y01i <- which.max(y01W)
  y01 <- uSub[y01i, ]
  
  
  
  
  
  png(
    filename = file.path("imgs3", sprintf("%03dobs1.png", seed)),
    width = 1024+512, height = 1024)
  par(mar = c(0,0,0,0), bg = 'lightgrey')
  plot(
    NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
    xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
    ylab=NA, xlab=NA)
  text(
    x = mean(range(uSub[,dm[1]])), y = 40,
    sprintf("t = %.2f", 0/100),
    cex = 3)
  points(uSub[, dm], pch=".")
  points(matrix(y0Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="green")
  dev.off()
  
  # for (i in seq(0, 100, by=1)) {
  # # for (i in c(0, deltaI)) {
  #   cat(i,",",sep="")
  #   png(
  #     filename = file.path("imgs3", sprintf("%03dframe%04d.png", seed, i)),
  #     width = 1024+512, height = 1024)
  #   par(mar = c(0,0,0,0), bg = 'lightgrey')
  #   plot(
  #     NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
  #     xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
  #     ylab=NA, xlab=NA)
  #   text(
  #     x = mean(range(uSub[,dm[1]])), y = 40,
  #     sprintf("t = %.2f", i/100),
  #     cex = 3)
  #   #points(uSub[-(1:i), dm], pch=".")
  #   points(uSub[-(1:i), dm], pch=".", col = cols[y0ColI])
  #   points(matrix(y0Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="green")
  #   dev.off()
  # }
  # 
  # png(
  #   filename = file.path("imgs3", sprintf("%03dobs2a.png", seed)),
  #   width = 1024+512, height = 1024)
  # par(mar = c(0,0,0,0), bg = 'lightgrey')
  # plot(
  #   NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
  #   xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
  #   ylab=NA, xlab=NA)
  # text(
  #   x = mean(range(uSub[,dm[1]])), y = 40,
  #   sprintf("t = %.2f", deltaI/100),
  #   cex = 3)
  # points(uSub[-(1:deltaI), dm], pch=".", col = cols[y0ColI])
  # points(matrix(y0Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="green")
  # points(matrix(y1Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="blue")
  # dev.off()
  # 
  # 
  # png(
  #   filename = file.path("imgs3", sprintf("%03dobs2b.png", seed)),
  #   width = 1024+512, height = 1024)
  # par(mar = c(0,0,0,0), bg = 'lightgrey')
  # plot(
  #   NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
  #   xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
  #   ylab=NA, xlab=NA)
  # text(
  #   x = mean(range(uSub[,dm[1]])), y = 40,
  #   sprintf("t = %.2f", deltaI/100),
  #   cex = 3)
  # points(uSub[-(1:deltaI), dm], pch=".", col = cols[y1ColI])
  # points(matrix(y0Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="green")
  # points(matrix(y1Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="blue")
  # dev.off()
  # 
  # png(
  #   filename = file.path("imgs3", sprintf("%03dobs2c.png", seed)),
  #   width = 1024+512, height = 1024)
  # par(mar = c(0,0,0,0), bg = 'lightgrey')
  # plot(
  #   NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
  #   xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
  #   ylab=NA, xlab=NA)
  # text(
  #   x = mean(range(uSub[,dm[1]])), y = 40,
  #   sprintf("t = %.2f", deltaI/100),
  #   cex = 3)
  # points(uSub[-(1:deltaI), dm], pch=".", col = cols[y01ColI])
  # points(matrix(y0Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="green")
  # points(matrix(y1Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="blue")
  # points(matrix(uSub[y01i+deltaI, dm], nrow=1), pch=4, lwd=5, cex=5, col="red")
  # dev.off()
  # 
  # png(
  #   filename = file.path("imgs3", sprintf("%03dobs2d.png", seed)),
  #   width = 1024+512, height = 1024)
  # par(mar = c(0,0,0,0), bg = 'lightgrey')
  # plot(
  #   NA, xaxt='n', yaxt='n', ann=FALSE, frame.plot=FALSE,
  #   xlim=range(uSub[,dm[1]]), ylim=range(uSub[,dm[2]]),
  #   ylab=NA, xlab=NA)
  # text(
  #   x = mean(range(uSub[,dm[1]])), y = 40,
  #   sprintf("t = %.2f", 0/100),
  #   cex = 3)
  # points(uSub[1:nSub, dm], pch=".", col = cols[y01ColI])
  # points(matrix(y0Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="green")
  # points(matrix(y1Proj[dm], nrow=1), pch=4, lwd=5, cex=5, col="blue")
  # points(matrix(y01[dm], nrow=1), pch=4, lwd=5, cex=5, col="red")
  # dev.off()

}