# express uncertainty in measurement of y1 by sampling "particles" around y1 with normal distri
# a particle filter
estimate <- function(attractor, y0, y1, timedNnFun, sd, particleCount) {
  cat("\t\t\t\t1")
  noise <- matrix(rnorm(particleCount*3, sd=sd), nrow=particleCount)
  cat("2")
  particles <- rep(y1, each = particleCount) + noise
  cat("3")
  idxes <- apply(particles, 1, \(p) timedNnFun(p)$idx)
  cat("4")
  tabulation <- tabulate(idxes)
  cat("5")
  nonZero <- tabulation != 0
  cat("6")
  estiIdx <- which(nonZero)
  cat("7")
  tabWeights <- tabulation[nonZero]
  cat("8")
  esti <- attractor$u[estiIdx, , drop=FALSE]
  cat("9")
  dstEsti <- sqrt(rowSums((esti - rep(y0, each = nrow(esti)))^2))
  cat("a")
  y0weights <- dnorm(dstEsti, sd=sd)
  cat("b")
  w <- tabWeights * y0weights
  cat("c")
  
  # mean estimate:
  # meanEsti <- colSums(esti * w) / sum(w)
  
  # TODO: make bandwidth (and evtl kernel a hyperparameter) or choose adaptively
  # TODO: kernel density like this is very memory demanding. Use KNN?
  # MAP estimate:
  dst <- as.matrix(stats::dist(esti))
  cat("d")
  kernel <- function(x) ifelse(x < 1, 1 - x^2, 0)
  cat("e")
  bandwidth <- 1
  cat("f")
  liklihood <- colSums(kernel(dst / bandwidth) * w)
  cat("g")
  map <- esti[which.max(liklihood), ]
  
  cat("h\n")
  
  return(map)
}

observeAndEstimate <- function(attractor, x0, x1, timedNnFun, sd, particleCount) {
  
  pt <- proc.time()
  cat("\t\t\tsample observations\n")
  # add noise to create observations
  y0 <- x0 + rnorm(3, sd=sd)
  y02 <- x0 + rnorm(3, sd=sd)
  y1 <- x1 + rnorm(3, sd=sd)

  cat("\t\t\testimate\n")
  # estimate
  estWait <- estimate(attractor, y0, y1, timedNnFun, sd, particleCount)
  estWaitPrj <- attractor$u[attractor$nnFun(estWait)$idx, , drop=FALSE]
  estNow <- (y0 + y02) / 2
  estNowPrj <- attractor$u[attractor$nnFun(estNow)$idx, , drop=FALSE]
  
  cat("\t\t\tduration: ", sprintf("%.1fs", (proc.time()-pt)[3]),"\n")
  
  # return squared distance to truth
  c(wait = sum((estWaitPrj-x0)^2),
    now = sum((estNowPrj -x0)^2))
}



generateAndError <- function(attractor, deltaI, timedNnFun, sd, noiseReps, particleCount) {
  
  cat("\t\tsample truth\n")
  # randomly draw the true locations of the two measurements
  i0 <- sample.int(attractor$n - deltaI, 1)
  i1 <- i0 + deltaI
  x0 <- attractor$u[i0, ]
  x1 <- attractor$u[i1, ]
  
  cat("\t\tobserveAndEstimate() with noiseReps = ", noiseReps, "\n")
  sqrErrors <- replicate(
    noiseReps, 
    observeAndEstimate(attractor, x0, x1, timedNnFun, sd, particleCount))
  
  return(sqrErrors)
}
  
run <- function(attractor, deltaT, sd, noiseReps, locationReps, particleCount) {
  
  deltaI <- round(deltaT / attractor$tStep)
  
  cat("\tbuild timedNnFun\n") # TODO: do not do this twice (see nnFun)
  timedNnFun <- FastKNN::buildKnnFunction(
    attractor$u[(deltaI+1):attractor$n,], 
    k = 1,
    removeNaRows = FALSE)
  
  cat("\tgenerateAndError() with locationReps = ", locationReps, "\n")
  sqrErrorss <- replicate(
    locationReps, 
    generateAndError(attractor, deltaI, timedNnFun, sd, noiseReps, particleCount))
  
  cat("\tclean up timedNnFun\n")
  FastKNN::deleteQueryFunction(timedNnFun)
  
  return(sqrErrorss)
}

execute <- function(
  deltaT = NULL,
  sd = NULL,
  noiseReps = NULL,
  locationReps = NULL,
  particleCount = NULL,
  outFile = NULL,
  attrFile = NULL
) {
  
  if (is.null(deltaT)) deltaT <- 0
  if (is.null(sd)) sd <- 1
  if (is.null(noiseReps)) noiseReps <- 1
  if (is.null(locationReps)) locationReps <- 100
  if (is.null(particleCount)) particleCount <- 1e4
  if (is.null(outFile)) outFile <- paste0("results_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RDS")
  if (is.null(attrFile)) attrFile <- "attractorLorenz63.RDS"
  
  cat("read attractor file\n")
  attractor <- readRDS(attrFile)
  
  cat("build nnFun\n")
  nnFun <- FastKNN::buildKnnFunction(
    attractor$u, 
    k = 1,
    removeNaRows = FALSE)
  attractor$nnFun <- nnFun
  
  cat("run\n")
  res <- run(
    attractor, 
    deltaT = deltaT, 
    sd = sd, 
    noiseReps = noiseReps, 
    locationReps = locationReps, 
    particleCount = particleCount)
  
  cat("save results\n")
  saveRDS(res, file = outFile)
  
  cat("clean up nnFun\n")
  FastKNN::deleteQueryFunction(nnFun)
}

args <- commandArgs(TRUE)
if (length(args) > 0) {
  argMat <- matrix(unlist(strsplit(args, "=")), nrow=2)
  numericValues <- as.numeric(argMat[2,])
  argList <- as.list(ifelse(is.na(numericValues), argMat[2,], numericValues))
  names(argList) <- argMat[1,]
} else {
  argList <- list()
}

execute(
  deltaT = argList$deltaT,
  sd = argList$sd,
  noiseReps = argList$noiseReps,
  locationReps = argList$locationReps,
  particleCount = argList$particleCount,
  outFile = argList$outFile,
  attrFile = argList$attrFile)
