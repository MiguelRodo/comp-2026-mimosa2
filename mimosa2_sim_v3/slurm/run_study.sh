#!/bin/bash
# =============================================================================
# Run ONE study on ONE node.     Usage (from the mimosa2_sim_v3 folder):
#   sbatch slurm/run_study.sh baseline        (or prior, heterogeneity, imbalance, maxit_check)
#   SIM_PROFILE=extended sbatch slurm/run_study.sh prior
# Re-submitting the same command after a time-out CONTINUES where it stopped.
# EDIT the #SBATCH lines and the module line to match your account / cluster
# (copy them from the job script you used for version 2).
# =============================================================================
#SBATCH --job-name=mimosa2_sim
#SBATCH --nodes=1
#SBATCH --ntasks=20                 # = number of parallel workers (one node)
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --output=slurm_%x_%j.out
##SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=ada

STUDY=${1:?"give the study name: baseline | prior | heterogeneity | imbalance | maxit_check"}
export SIM_PROFILE=${SIM_PROFILE:-standard}
export N_WORKERS=${SLURM_NTASKS:-20}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# [CHANGE 5 Oct] Your cluster runs R inside the "mimosa2" Apptainer container
# (as in your sim.sh), not the system R. run_r runs one R script in it.
run_r() { apptainer-rscript -f mimosa2 -- "source(\"$1\")"; }
export R_LIBS=/scratch/abrmoe030/R_libs${R_LIBS:+:$R_LIBS}   # your package library (all R processes)        # <- same module line as before

case $STUDY in
  baseline)      SCRIPT=sims/20_sim_baseline.R ;;
  prior)         SCRIPT=sims/21_sim_prior.R ;;
  heterogeneity) SCRIPT=sims/22_sim_heterogeneity.R ;;
  imbalance)     SCRIPT=sims/23_sim_imbalance.R ;;
  maxit_check)   SCRIPT=sims/24_sim_maxit_check.R ;;
  *) echo "unknown study $STUDY"; exit 1 ;;
esac
echo "$(date) start $STUDY profile=$SIM_PROFILE workers=$N_WORKERS"
run_r "$SCRIPT"
echo "$(date) end $STUDY"
