#!/usr/bin/env bash
# =============================================================================
# OPTION 2 (more cores): split one study over several nodes with a job array.
# Each array element is a separate 20-core (or larger) job that runs every
# K-th task; all write to the same _simulations/<study>/tasks folder.
#   SIM_PROFILE=extended sbatch --array=0-3 slurm/run_study_array.sh baseline
#   ... repeat for prior, heterogeneity, imbalance, maxit_check ...
# When ALL array jobs of a study have finished, combine them:
#   SIM_PROFILE=extended sbatch --dependency=afterany:<ARRAY_JOB_ID> slurm/combine.sh baseline
# (or run  SIM_PROFILE=extended Rscript sims/29_combine.R baseline  on a login/interactive node).
# Array IDs must start at 0, or SLURM_ARRAY_TASK_MIN is subtracted automatically.
# Re-submitting the same array after a time-out continues where it stopped.
# =============================================================================
#SBATCH --job-name=mimosa2_arr
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --output=slurm_%x_%A_%a.out
##SBATCH --account=YOUR_ACCOUNT
##SBATCH --partition=YOUR_PARTITION

STUDY=${1:?"give the study name"}
export SIM_PROFILE=${SIM_PROFILE:-extended}
export N_WORKERS=${SLURM_NTASKS:-20}
# number of chunks = number of array elements (older SLURM versions do not set
# SLURM_ARRAY_TASK_COUNT: then submit with  N_CHUNKS=4 sbatch --array=0-3 ...)
export N_CHUNKS=${N_CHUNKS:-${SLURM_ARRAY_TASK_COUNT:?"submit with --array, or set N_CHUNKS"}}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# module load software/R-4.x.x

case $STUDY in
  baseline)      SCRIPT=sims/20_sim_baseline.R ;;
  prior)         SCRIPT=sims/21_sim_prior.R ;;
  heterogeneity) SCRIPT=sims/22_sim_heterogeneity.R ;;
  imbalance)     SCRIPT=sims/23_sim_imbalance.R ;;
  maxit_check)   SCRIPT=sims/24_sim_maxit_check.R ;;
  *) echo "unknown study $STUDY"; exit 1 ;;
esac
echo "$(date) start $STUDY chunk $SLURM_ARRAY_TASK_ID of $N_CHUNKS"
Rscript "$SCRIPT"
echo "$(date) end"
