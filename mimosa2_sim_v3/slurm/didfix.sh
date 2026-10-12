#!/bin/bash
# =============================================================================
# didfix.sh : apply the v3.1 DiD continuity correction to the EXISTING results
# files, for both profiles. No MIMOSA2 fit is repeated. No tables/figures are
# made (run the analysis scripts afterwards, with the profile each study uses).
#
# Usage (from the mimosa2_sim_v3 folder):
#   sbatch slurm/didfix.sh
#
# Needs these files in _simulations/<study>/ :
#   extended: heterogeneity, imbalance, prior   (<study>_results_extended.rds)
#   standard: baseline, maxit_check (+ any other standard files present)
# Missing files are skipped with "not found"; already-corrected files are
# skipped with "already corrected". Each corrected file gets a backup
# <study>_results_<profile>_before_didfix.rds first.
# =============================================================================
#SBATCH --job-name=mimosa2_didfix
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --output=slurm_%x_%j.out
##SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=ada

export N_WORKERS=${SLURM_NTASKS:-20}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run_r() { apptainer-rscript -f mimosa2 -- "source(\"$1\")"; }
export R_LIBS=/scratch/abrmoe030/R_libs${R_LIBS:+:$R_LIBS}

set -e
cd /scratch/abrmoe030/projects/mimosa2/mimosa2_sim_v3

echo "$(date) results files found:"
ls -la _simulations/*/*_results_*.rds

# ---- 1. extended profile: the three studies run with the extended grid ----
echo "$(date) ==== DiD fix: extended ===="
export SIM_PROFILE=extended
export STUDIES="heterogeneity imbalance prior"
run_r sims/40_recompute_did.R

# ---- 2. standard profile: baseline and maxit_check (and the other standard
#         files too, so every results file on disk is corrected) ----
echo "$(date) ==== DiD fix: standard ===="
export SIM_PROFILE=standard
export STUDIES="baseline maxit_check heterogeneity imbalance prior"
run_r sims/40_recompute_did.R

echo "$(date) done"
