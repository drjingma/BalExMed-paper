#!/bin/bash
## =====================================================================
## SLURM array for the simulation study of Section 4.
## One array task runs a chunk of 20 replicates of one grid cell in parallel
## across cpus-per-task cores: 42 cells x ceil(200 / 20) = 420 tasks.
## Each replicate is written to its own file and skipped if present, so a
## cancelled or timed-out task loses nothing: resubmit and finished replicates
## are skipped.
##
##   submit (from simulations/):  sbatch submit_sim.sh
##   collect:                     Rscript sim_run.R aggregate
##   tables and figures:          Rscript sim_report.R
## If CHUNK changes, reset --array from:  Rscript sim_run.R ntasks <CHUNK>
## Adapt the module line and resource requests to your cluster.
## =====================================================================

#SBATCH --job-name=balexmed_sim
#SBATCH --array=0-419            # == Rscript sim_run.R ntasks 20
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=1-0
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err

set -euo pipefail
# module load R                  # load R as required on your cluster
mkdir -p logs
CHUNK=20                         # must match the chunk used to size --array
cd "${SLURM_SUBMIT_DIR}"
echo "host=$(hostname) task=$SLURM_ARRAY_TASK_ID cpus=$SLURM_CPUS_PER_TASK chunk=$CHUNK"
Rscript sim_run.R task "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK" "$CHUNK"
