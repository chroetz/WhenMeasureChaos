files <- dir("./data3", pattern = "^WhenMeasureChaos_\\d+_.*\\.RDS")
paths <- paste0("./data3/", files)

n <- length(files)
sds <- double(n)
deltaTs <- double(n)

attractor <- readRDS("C:/Users/cschoetz/Documents/WhenMeasureChaos/attractorLorenz63.RDS")

mses <- function(idxes) {
  apply(idxes, 3, \(i) {
    x0 <- attractor$u[i[1,1],]
    esti1 <- attractor$u[i[1,-1],,drop=FALSE]
    esti2 <- attractor$u[i[2,-1],,drop=FALSE]
    se1 <- rowSums((rep(x0, each=nrow(esti1)) - esti1)^2)
    se2 <- rowSums((rep(x0, each=nrow(esti2)) - esti2)^2)
    c(mean(se1), mean(se2), sd(se1), sd(se2))})
}

# TODO: unclear..
trajLenHalf <- 1e3
projErr <- function(idxes) {
  apply(idxes, 3, \(i) {
    i0 <- i[1,1]
    x0 <- attractor$u[i0,]
    trajI <- pmin(attractor$n, pmax(1, (i0 - trajLenHalf):(i0 + trajLenHalf)))
    traj <- attractor$u[trajI, ]
    
    # TODO: this is not a fair projection as traj could be a whole circle
    dst <- rowSums((rep(x0, each=nrow(traj)) - traj)^2)
    ddst <- diff(dst)
    left <- which.max(c(rev(ddst[1:trajLenHalf] > 0), TRUE))
    right <- which.max(c(ddst[(trajLenHalf+1):(2*trajLenHalf)] < 0, TRUE))
    len <- min(c(right, left))
    leftTraj <-  traj[trajLenHalf:(trajLenHalf+1-len),]
    rightTraj <-  traj[(trajLenHalf+2):(trajLenHalf+1+len),]
    dst <- rowSums((leftTraj-rightTraj)^2)
    llen <- which.max(c(diff(dst) < 0, TRUE))
    traj <- traj[(trajLenHalf+1-llen):(trajLenHalf+1+llen),]
    trajNnFun <- FastKNN::buildKnnFunction(traj, 1)
    
    esti1 <- attractor$u[i[1,-1],,drop=FALSE] # MAP
    #esti2 <- attractor$u[i[2,-1],,drop=FALSE] # single
    
    dsts <- sapply(1:nrow(esti1), \(j) {
      nn <- trajNnFun(esti1[j,])
      c(abs(left+1 - nn$idx), nn$distSqr)
    })
    return(dsts)
  })
}

ses <- function(idxes) {
  apply(idxes, 3, \(i) {
    x0 <- attractor$u[i[1,1],]
    esti1 <- attractor$u[i[1,-1],,drop=FALSE]
    esti2 <- attractor$u[i[2,-1],,drop=FALSE]
    se1 <- rowSums((rep(x0, each=nrow(esti1)) - esti1)^2)
    se2 <- rowSums((rep(x0, each=nrow(esti2)) - esti2)^2)
    cbind(se1, se2)})
}

# sess <- sapply(seq_len(n), \(i) {
#   rds <- readRDS(paths[i])
#   ses(rds$idxes)
#   })
# dim(sess) <- c(100, 2, 1000, 36)

msesLst <- list()
projDstLst <- list()
for (i in seq_len(n)) {
  cat(i, ",", sep="")
  rds <- readRDS(paths[i])
  sds[i] <- rds$opts$sd
  deltaTs[i] <- rds$opts$deltaT
  msesLst[[i]] <- mses(rds$idxes)
  projDstLst[[i]] <- projErr(rds$idxes)
}

projDstArray <- unlist(projDstLst)
dim(projDstArray) <- c(2, 100, 200, length(projDstLst))

geomMean <- sapply(msesLst, \(x) exp(mean(log(x[2,]/x[1,]))))
geomVar <- sapply(msesLst, \(x) mean((exp(log(x[2,]/x[1,]) - mean(log(x[2,]/x[1,]))))^2))
arithMean <- sapply(msesLst, \(x) mean(x[1,]))
mmse1 <- sapply(msesLst, \(x) mean(x[2,]))
arithVar <- sapply(msesLst, \(x) var(x[1,]))
meanVar <- sapply(msesLst, \(x) mean(x[3,]^2))

#totalVar <- apply(sess, c(2,4), \(x) var(as.vector(x)))
totalVar <- meanVar + arithVar
  
library(tidyverse)
data <- 
  tibble(
    sd = sds,
    deltaT = deltaTs,
    improvement = geomMean,
    uncertainty = geomVar,
    meanMSE = arithMean,
    locationVar = arithVar,
    noiseVariance = meanVar,
    totalVar = totalVar,
    mmse1 = mmse1,
    attrDst = apply(projDstArray[2,,,], 3, mean),
    timeDst = apply(projDstArray[1,,,], 3, mean)
  )
data <- 
  data |> 
  mutate(sd = factor(sd))
data |> 
  ggplot() + geom_line(aes(x = deltaT, y = improvement, color = sd))
pltMse <- 
  data |>
  group_by(sd) |> 
  arrange(deltaT) |> 
  mutate(meanMSE = meanMSE/meanMSE[1]) |> 
  ggplot() + 
  geom_line(aes(x = deltaT, y = meanMSE, color = sd)) +
  labs(x="t", y="MSE", color=expression(sigma)) +
  geom_abline(slope = 0, intercept=1, color = "black") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0.6, 2.2, by = 0.2))
pltMsePerp <- 
  data |>
  group_by(sd) |> 
  arrange(deltaT) |> 
  mutate(attrDst = attrDst/attrDst[1]) |> 
  ggplot() + 
  geom_line(aes(x = deltaT, y = attrDst, color = sd)) +
  labs(x="t", y="perpendicular MSE", color=expression(sigma)) +
  geom_abline(slope = 0, intercept=1, color = "black") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0.0, 2.2, by = 0.2))
pltMsePara <- 
  data |>
  group_by(sd) |> 
  arrange(deltaT) |> 
  mutate(timeDst = timeDst/timeDst[1]) |> 
  ggplot() + 
  geom_line(aes(x = deltaT, y = timeDst, color = sd)) +
  labs(x="t", y="parallel MSE", color=expression(sigma)) +
  geom_abline(slope = 0, intercept=1, color = "black") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0.0, 2.2, by = 0.1))
data |> 
  filter(sd == 0.6) |> 
  mutate(
    low = meanMSE - 2*sqrt(totalVar / 2e4),
    high = meanMSE + 2*sqrt(totalVar / 2e4)
  ) |> 
  ggplot(aes(x = deltaT, y = meanMSE, ymin = low, ymax = high)) +
  geom_ribbon(alpha=0.3) +
  geom_line() + 
  geom_line(aes(x = deltaT, y = mmse1/2), color = "red") + ylim(c(0,NA))
