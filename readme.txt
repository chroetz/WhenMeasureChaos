cd into this folder
Rscript ./startSlurm.R # starts run.sh which executes module load ./Rmodule and Rscript lorenz63run.R
squeue -u cschoetz # check jobs
scancel <jobid> # cancel jobs
