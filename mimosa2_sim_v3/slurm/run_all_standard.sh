#!/bin/bash
# =============================================================================
# OPTION 1 (about 20 cores): all studies one after another in ONE job.
#   sbatch slurm/run_all_standard.sh
# Order: the four smaller studies first (about 8 h on 20 cores), then the long
# baseline (about 11 h). The baseline runs repetition by repetition, so if
# time runs out you can stop it and still analyse every scenario with the
# repetitions completed so far (sims/29_combine.R uses complete repetitions). Every study is resumable, so if
# the job hits its time limit just submit it again: finished tasks are skipped.
# Projected time is printed by sims/11_smoke_test.R (section 5). If your
# cluster lets you run several jobs at once, submitting the studies
# separately with slurm/run_study.sh is faster.
# =============================================================================
#SBATCH --job-name=mimosa2_all
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --output=slurm_%x_%j.out
##SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=ada

export SIM_PROFILE=standard
export N_WORKERS=${SLURM_NTASKS:-20}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# [CHANGE 5 Oct] Your cluster runs R inside the "mimosa2" Apptainer container
# (as in your sim.sh), not the system R. run_r runs one R script in it.
run_r() { apptainer-rscript -f mimosa2 -- "source(\"$1\")"; }
export R_LIBS=/scratch/abrmoe030/R_libs${R_LIBS:+:$R_LIBS}   # your package library (all R processes)

set -e
run_r sims/10_dgm_checks.R
for S in sims/21_sim_prior.R sims/23_sim_imbalance.R sims/22_sim_heterogeneity.R sims/24_sim_maxit_check.R sims/20_sim_baseline.R; do
  echo "$(date) start $S"; run_r "$S"; echo "$(date) end $S"
done
run_r analysis/30_performance_tables.R
