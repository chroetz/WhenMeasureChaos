args <- commandArgs(TRUE)
argMat <- trimws(matrix(unlist(strsplit(args, "=")), nrow=2))
argList <- strsplit(argMat[2,], ",")
grid <- expand.grid(argList, stringsAsFactors=FALSE)
# grid <- lapply(grid, \(column) if (!any(is.na(as.numeric(column)))) as.numeric(column) else column)
# grid <- as.data.frame(grid)
names(grid) <- argMat[1,]
pasteRow <- function(row) {
  paste0(names(row), "=", as.character(row), collapse=" ")
}

scriptName <- "run.R"

for (i in seq_len(nrow(grid))) {
  jobName <- paste0("WhenMeasureChaos_", i, "_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
  cat("Starting SLURM job", jobName, "\n")
  clcom <- paste0(
    "sbatch",
    " --qos=short",
    " --ntasks=1",
    " --cpus-per-task=4", # RAM: ~ cpus-per-task * 4 GB
    " --job-name=", jobName,
    " --output=", jobName, "_%j.out",
    " --error=", jobName, "_%j.err",
    " --mail-type=END",
    " --wrap=\"Rscript ", scriptName, " ", pasteRow(grid[i,]), "\"")
  cat(clcom, "\n")
  system(clcom)
}
