#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name="sim"
#SBATCH --partition=ada
start_time=$(date +%s)

echo "HOSTNAME: $HOSTNAME"

echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "SLURM_NTASKS=$SLURM_NTASKS"
echo "SLURM_JOB_CPUS_PER_NODE=$SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE"

echo "SLURM CPU affinity:"
taskset -pc $$

echo "Available processors:"
nproc

echo " "
echo " "
echo " "
echo "-------------------"
date

apptainer-rscript -f mimosa2 -- 'source("_tmp/Non_beta_simulations.R"); source("_tmp/Prior_simulations.R")'
echo "Completed running simulation"
date
echo "-------------------"
echo " "

# Record the end time
end_time=$(date +%s)

# Calculate the duration
duration=$((end_time - start_time))

# Convert duration to human-readable format
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

# Append the duration to the Slurm standard output log
echo "--- Script Duration ---"
printf "Elapsed time: %02d:%02d:%02d\n" $hours $minutes $seconds
echo "-----------------------"
