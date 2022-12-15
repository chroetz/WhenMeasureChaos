files <- dir("./data5", pattern = "^WhenMeasureChaos_\\d+_.*\\.RDS")
paths <- paste0("./data5/", files)

# paths <- paths[1]

n <- length(paths)
sds <- double(n)
deltaTs <- double(n)

attractor <- readRDS("C:/Users/cschoetz/Documents/WhenMeasureChaos/attractorLorenz63.RDS")

mses <- function(idxes) {
  x0 <- attractor$u[idxes[1,],]
  idxes[2, idxes[2,] == 0] <- NA_integer_
  esti <- attractor$u[idxes[2,],]
  dstSq <- rowSums((x0 - esti)^2)
  dstSq
}

trajLenHalf <- 1e3
projErr <- function(idxes) {
  apply(idxes, 2, \(i) {
    
    # The truth.
    i0 <- i[1]
    x0 <- attractor$u[i0,]
    
    # Create a trajectory segment around the truth such that:
    # * When moving along the segment in one direction away from the truth, (space-)distance to the truth does not decrease.
    # * When moving along the segment in both direction away from the truth with the same speed, (space-)distance between the end-points increases.
    
    trajI <- pmin(attractor$n, pmax(1, (i0 - trajLenHalf):(i0 + trajLenHalf)))
    traj <- attractor$u[trajI, ]
    
    dst <- DEEBesti:::distSqrToVec(traj, x0)
    ddst <- diff(dst)
    left <- which.max(c(rev(ddst[1:trajLenHalf] > .Machine$double.eps), TRUE))
    right <- which.max(c(ddst[(trajLenHalf+1):(2*trajLenHalf)] < -.Machine$double.eps, TRUE))
    len <- min(c(right, left))
    stopifnot(len > 1)
    leftTraj <- traj[trajLenHalf:(trajLenHalf+1-len),]
    rightTraj <- traj[(trajLenHalf+2):(trajLenHalf+1+len),]
    dst <- rowSums((leftTraj-rightTraj)^2)
    llen <- which.max(c(diff(dst) < 0, TRUE))
    trajShort <- traj[(trajLenHalf+1-llen):(trajLenHalf+1+llen),]
    # trajShort is the final trajectory segment
    
    # Calculate error measures.
    esti1 <- attractor$u[i[2],] # MAP
    if (length(esti1) == 0) return(rbind(NA_real_,NA_real_,NA_real_))
    dstSqEsti <- DEEBesti:::distSqrToVec(trajShort, esti1)
    nnI <- which.min(dstSqEsti)
    truthI <- llen+1
    nnDstSq <- dstSqEsti[nnI]
    trajDstSq <- sum((trajShort[nnI, ] - trajShort[truthI,])^2)
    dsts <- rbind(abs(truthI - nnI), trajDstSq, nnDstSq)
    
    return(dsts)
  })
}

msesLst <- list()
projDstLst <- list()

for (i in (length(msesLst)+1):n) {
  cat(i, ",", sep="")
  rds <- readRDS(paths[i])
  sds[i] <- rds$opts$sd
  deltaTs[i] <- rds$opts$deltaT
  msesLst[[i]] <- mses(rds$idxes)
  projDstLst[[i]] <- projErr(rds$idxes)
  
  # saveRDS(list(sds = sds, deltaTs= deltaTs, msesLst=msesLst, projDstLst=projDstLst),
  #         file=sprintf("dump%03d.RDS", i))
}
saveRDS(
  list(
    sds = sds, 
    deltaTs = deltaTs, 
    msesLst = msesLst, 
    projDstLst = projDstLst),
  file=sprintf("data5-%d.RDS", Sys.time() |> as.integer()))
