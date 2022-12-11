estimate <- function(attractor, y0, y1, nnFun, sd, deltaI) {
  i0 <- nnFun(y0)$idx
  i1 <- i0 + deltaI
  esti0 <- attractor$u[i0, , drop=FALSE]
  esti1 <- attractor$u[i1, , drop=FALSE]
  dstEsti0 <- sqrt(rowSums((esti0 - rep(y0, each = nrow(esti0)))^2))
  dstEsti1 <- sqrt(rowSums((esti1 - rep(y1, each = nrow(esti1)))^2))
  w0 <- dnorm(dstEsti0, sd=sd)
  w1 <- dnorm(dstEsti1, sd=sd)
  mapI <- i0[which.max(w0 * w1)]
  return(c(mapI, i0[1]))
}


observeAndEstimate <- function(
    attractor, x0, x1, nnFun, sd, deltaI
) {
  
  #pt <- proc.time()
  logg(3, "sample observations")
  # add noise to create observations
  y0 <- x0 + rnorm(3, sd=sd)
  y1 <- x1 + rnorm(3, sd=sd)

  logg(3, "estimate")
  # estimate
  estiI <- estimate(attractor, y0, y1, nnFun, sd, deltaI)
  #logg(3, "duration: ", sprintf("%.1fs", (proc.time()-pt)[3]))
  
  # return estimated index
  estiI
}


generateAndError <- function(
    attractor, deltaI, nnFun, sd, noiseReps
) {
  
  logg(2, "sample truth")
  # randomly draw the true locations of the two measurements
  i0 <- sample.int(attractor$n - deltaI, 1)
  i1 <- i0 + deltaI
  x0 <- attractor$u[i0, ]
  x1 <- attractor$u[i1, ]
  
  logg(2, "observeAndEstimate() with noiseReps = ", noiseReps)
  estiIs <- replicate(
    noiseReps, 
    observeAndEstimate(attractor, x0, x1, nnFun, sd, deltaI))
  
  return(cbind(truth = i0, estiIs))
}

  
run <- function(
    attractor, deltaT, sd, noiseReps, locationReps, particleCount
) {
  
  deltaI <- round(deltaT / attractor$tStep)
  
  logg(1, "build nnFun")
  nnFun <- FastKNN::buildKnnFunction(
    attractor$u[1:(attractor$n - deltaI), ], 
    k = particleCount,
    removeNaRows = FALSE)
  
  logg(1, "generateAndError() with locationReps = ", locationReps)
  idxes <- replicate(
    locationReps, 
    generateAndError(attractor, deltaI, nnFun, sd, noiseReps))
  
  logg(1, "clean up nnFun")
  FastKNN::deleteQueryFunction(nnFun)
  
  return(idxes)
}


execute <- function(opts) {
  
  # Set default opts.
  if (is.null(opts$deltaT)) opts$deltaT <- 0
  if (is.null(opts$sd)) opts$sd <- 1
  if (is.null(opts$noiseReps)) opts$noiseReps <- 1
  if (is.null(opts$locationReps)) opts$locationReps <- 100
  if (is.null(opts$particleCount)) opts$particleCount <- 1e4
  if (is.null(opts$outFile)) opts$outFile <- paste0("results_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".RDS")
  if (is.null(opts$attrFile)) opts$attrFile <- "attractorLorenz63.RDS"
  
  logg(0, "Use following opts:")
  for (i in seq_along(opts)) 
    logg(0, "* ", names(opts)[i], ": ", opts[[i]], " (", typeof(opts[[i]]),")\n", sep="")

  logg(0, "read attractor file")
  attractor <- readRDS(opts$attrFile)
  
  logg(0, "run")
  idxes <- run(
    attractor, 
    deltaT = opts$deltaT, 
    sd = opts$sd, 
    noiseReps = opts$noiseReps, 
    locationReps = opts$locationReps, 
    particleCount = opts$particleCount)
  
  logg(0, "save results")
  res <- list(
    idxes = idxes,
    opts = opts)
  saveRDS(
    res, 
    file = opts$outFile)
  return(res)
}


logg <- function(level, ...) {
  if (level <= .logThreshold) {
    cat(paste0(c(rep("\t", level), ..., "\n"), collapse=""))
  }
}

.logThreshold <- 2



args <- commandArgs(TRUE)
if (interactive()) {
  args <- c(
    "particleCount = 1e5",
    "noiseReps = 20",
    "locationReps = 1",
    "sd = 0.1",
    "deltaT = 3")
}

if (length(args) > 0) {
  argMat <- matrix(trimws(unlist(strsplit(args, "="))), nrow=2)
  numericValues <- suppressWarnings(as.numeric(argMat[2,]))
  argList <- ifelse(is.na(numericValues), as.list(argMat[2,]), as.list(numericValues))
  names(argList) <- argMat[1,]
} else {
  argList <- list()
}

pt <- proc.time()
prof <- profvis::profvis(
res <- execute(argList))
print(proc.time() - pt)
