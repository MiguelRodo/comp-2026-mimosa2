#!/bin/bash
# =============================================================================
# v3.1 (6 Oct 2026): apply the DiD fix to the EXISTING results and redo every
# table and figure. No MIMOSA2 fit is repeated.  Usage (mimosa2_sim_v3 folder):
#   sbatch slurm/fix_v31.sh
# If the baseline is being re-run with the v3.1 code, submit this job to start
# after it, so the baseline figures are made too:
#   JOB=$(sbatch --parsable slurm/run_study.sh baseline)
#   sbatch --dependency=afterok:$JOB slurm/fix_v31.sh
# A study run with the v3.1 code already has the corrected DiD and is skipped
# by 40_recompute_did.R automatically.
# =============================================================================
#SBATCH --job-name=mimosa2_fix31
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=03:00:00
#SBATCH --mem=48G
#SBATCH --output=slurm_%x_%j.out
##SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=ada

export SIM_PROFILE=standard
export N_WORKERS=${SLURM_NTASKS:-20}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run_r() { apptainer-rscript -f mimosa2 -- "source(\"$1\")"; }
export R_LIBS=/scratch/abrmoe030/R_libs${R_LIBS:+:$R_LIBS}

set -e
ls -la _simulations/*/*_results_standard.rds
run_r sims/40_recompute_did.R                 # all studies; skips ones already corrected; stops if a dataset does not reproduce
for S in analysis/30_performance_tables.R analysis/31_plots_baseline.R analysis/32_plots_prior.R \
         analysis/33_plots_heterogeneity.R analysis/34_plots_imbalance.R analysis/35_maxit_check.R; do
  echo "$(date) $S"; run_r "$S"
done
echo "$(date) done"
