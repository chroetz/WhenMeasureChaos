args <- commandArgs(TRUE)
#scriptName <- "test.R"
scriptName <- "lorenz63run.R"
jobName <- paste0("WhenMeasureChaos_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("Starting SLURM job", jobName, "\n")
clcom <- paste0(
  "sbatch",
  " --qos=short",
  " --ntasks=1",
  " --cpus-per-task=3", # RAM: ~ cpus-per-task * 4 GB
  " --job-name=", jobName,
  " --output=", jobName, "_%j.out",
  " --error=", jobName, "_%j.err",
  " --mail-type=END",
  " --wrap=\"Rscript ", scriptName, " ", paste(args, collapse=" "), "\"")
cat(clcom, "\n")
system(clcom)
