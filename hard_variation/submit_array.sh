#!/bin/bash
#SBATCH --job-name=thesis_array
#SBATCH --array=1-20
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/run_%A_%a.out
#SBATCH --error=logs/run_%A_%a.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

# one core per task: stop BLAS/OpenMP from spawning extra threads
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hsmm

Rscript run_array.R
