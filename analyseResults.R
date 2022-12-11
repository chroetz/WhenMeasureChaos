rds1 <- readRDS("data/results_2022-12-03_16-40-06.RDS")
files <- dir("./data", pattern = "^WhenMeasureChaos_\\d+_.*\\.RDS")
sel <- rep(FALSE, length(files))
sel <- sapply(files, \(f) {
  rds2 <- readRDS(paste0("data/", f))
  rds1$opts$sd == rds2$opts$sd && rds1$opts$deltaT == rds2$opts$deltaT
})
rds2 <- readRDS(paste0("data/", files[sel][2]))

attractor <- readRDS("attractorLorenz63.RDS")

mses <- function(idxes) {
  apply(idxes, 3, \(i) {
    x0 <- attractor$u[i[1,1],]
    esti1 <- attractor$u[i[1,-1],,drop=FALSE]
    esti2 <- attractor$u[i[2,-1],,drop=FALSE]
    se1 <- rowSums((rep(x0, each=nrow(esti1)) - esti1)^2)
    se2 <- rowSums((rep(x0, each=nrow(esti2)) - esti2)^2)
    c(mean(se1), mean(se2), sd(se1), sd(se2))})
}

mse1 <- mses(rds1$idxes)
mse2 <- mses(rds2$idxes)
mmse1 <- rowMeans(mse1)
mmse2 <- rowMeans(mse2)
