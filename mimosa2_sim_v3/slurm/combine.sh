#!/bin/bash
#SBATCH --job-name=mimosa2_combine
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=slurm_%x_%j.out
# Usage: SIM_PROFILE=extended sbatch slurm/combine.sh baseline prior heterogeneity imbalance maxit_check
export SIM_PROFILE=${SIM_PROFILE:-extended}
# [CHANGE 5 Oct] Your cluster runs R inside the "mimosa2" Apptainer container
# (as in your sim.sh), not the system R. run_r runs one R script in it.
run_r() { apptainer-rscript -f mimosa2 -- "source(\"$1\")"; }
export R_LIBS=/scratch/abrmoe030/R_libs${R_LIBS:+:$R_LIBS}   # your package library (all R processes)
export STUDIES="$*"
run_r sims/29_combine.R
run_r analysis/30_performance_tables.R
