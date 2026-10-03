#!/bin/bash
#SBATCH --job-name=mimosa2_combine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=slurm_%x_%j.out
# Usage: SIM_PROFILE=extended sbatch slurm/combine.sh baseline prior heterogeneity imbalance maxit_check
export SIM_PROFILE=${SIM_PROFILE:-extended}
# module load software/R-4.x.x
Rscript sims/29_combine.R "$@"
Rscript analysis/30_performance_tables.R
