args <- commandArgs(TRUE)
if (interactive()) {
  args <- c(
    "noiseReps = 1",
    "locationReps = 100",
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
res <- WhenMeasureChaosR::execute(argList) # 0.1s per rep
print(proc.time() - pt)
