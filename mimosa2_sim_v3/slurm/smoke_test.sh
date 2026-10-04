#!/bin/bash
# Run FIRST (about 15-30 minutes):   sbatch slurm/smoke_test.sh
#SBATCH --job-name=mimosa2_smoke
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --output=slurm_%x_%j.out
##SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=ada
export N_WORKERS=${SLURM_NTASKS:-8}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# [CHANGE 5 Oct] Your cluster runs R inside the "mimosa2" Apptainer container
# (as in your sim.sh), not the system R. run_r runs one R script in it.
run_r() { apptainer-rscript -f mimosa2 -- "source(\"$1\")"; }
export R_LIBS=/scratch/abrmoe030/R_libs${R_LIBS:+:$R_LIBS}   # your package library (all R processes)
run_r sims/10_dgm_checks.R && run_r sims/11_smoke_test.R
