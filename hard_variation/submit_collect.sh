#!/bin/bash
#SBATCH --job-name=thesis_collect
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/collect_%j.out
#SBATCH --error=logs/collect_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate hsmm

Rscript collect_runs.R
