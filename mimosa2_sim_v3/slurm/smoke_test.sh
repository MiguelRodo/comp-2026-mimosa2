#!/bin/bash
# Run FIRST (about 15-30 minutes):   sbatch slurm/smoke_test.sh
#SBATCH --job-name=mimosa2_smoke
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --output=slurm_%x_%j.out
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --partition=YOUR_PARTITION
export N_WORKERS=${SLURM_NTASKS:-8}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# module load software/R-4.x.x
Rscript sims/10_dgm_checks.R && Rscript sims/11_smoke_test.R
