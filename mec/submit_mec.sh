#!/bin/bash
## =====================================================================
## SLURM array for the eight MEC-APS fits in Figure 6: 8 settings x 4 chains,
## one single-core chain per array task. Run from mec/.
## One change at a time around the primary model (eta = 1/20, diabetes-adjusted,
## metformin users retained, total body fat as the adiposity covariate):
##   P   eta = 1/20, diabetes                        primary
##   A1  eta = 1/3,  diabetes                        prior
##   A2  eta = 1/5,  diabetes                        prior
##   A3  eta = 1/20, metformin                       adjustment
##   A4  eta = 1/20, neither                         adjustment
##   S2  eta = 1/20, diabetes, no adiposity covariate
##   S3  eta = 1/20, diabetes, metformin users excluded
##   S4  eta = 1/20, neither,  participants with diabetes excluded
##
##   submit all chains:   sbatch submit_mec.sh
##   combine chains:      bash submit_mec.sh aggregate   (login or interactive node)
##   rerun one task:      sbatch --array=<id> submit_mec.sh   (finished chains are skipped)
## Task id = 4 * (setting index, starting at 0) + (chain - 1).
## Adapt the module line and resource requests to your cluster.
## =====================================================================

#SBATCH --job-name=balexmed_mec
#SBATCH --array=0-31             # 8 settings x 4 chains
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --output=logs/mec_%A_%a.out
#SBATCH --error=logs/mec_%A_%a.err

set -euo pipefail
# module load R                  # load R as required on your cluster
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs

NCHAINS=4
## each setting: "eta eth adjust n_chains n_iter burn_in sampler covars sample"
SETTINGS=(
  "20 all diabet    4 5000 1000 gibbs new    retained"   # P
  "3  all diabet    4 5000 1000 gibbs new    retained"   # A1
  "5  all diabet    4 5000 1000 gibbs new    retained"   # A2
  "20 all metformin 4 5000 1000 gibbs new    retained"   # A3
  "20 all none      4 5000 1000 gibbs new    retained"   # A4
  "20 all diabet    4 5000 1000 gibbs noadip retained"   # S2
  "20 all diabet    4 5000 1000 gibbs new    nometf"     # S3
  "20 all none      4 5000 1000 gibbs new    not2d"      # S4
)
NSET=${#SETTINGS[@]}

if [ "${1:-}" = "aggregate" ]; then
  for cfg in "${SETTINGS[@]}"; do
    read ETA ETH ADJ NC NI BI SAM COV SMP <<< "$cfg"
    echo "== aggregate: $cfg =="
    Rscript 03_fit_mec.R "$ETA" "$ETH" "$ADJ" "$NC" "$NI" "$BI" "$SAM" aggregate "" "$COV" "$SMP"
  done
  exit 0
fi

S=$(( SLURM_ARRAY_TASK_ID / NCHAINS ))
K=$(( SLURM_ARRAY_TASK_ID % NCHAINS + 1 ))
if [ "$S" -ge "$NSET" ]; then echo "task $SLURM_ARRAY_TASK_ID beyond the settings list"; exit 0; fi
read ETA ETH ADJ NC NI BI SAM COV SMP <<< "${SETTINGS[$S]}"
echo "task $SLURM_ARRAY_TASK_ID -> setting $S {${SETTINGS[$S]}} chain $K"
Rscript 03_fit_mec.R "$ETA" "$ETH" "$ADJ" "$NC" "$NI" "$BI" "$SAM" chain "$K" "$COV" "$SMP"
