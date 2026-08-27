#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --job-name="sim"
#SBATCH --partition=ada

# Record the start time
start_time=$(date +%s)

echo "HOSTNAME: $HOSTNAME"

echo " "
echo " "
echo " "

echo "-------------------"
echo "Run simulation"
date
apptainer-rscript -f mimosa2 -- 'source("_tmp/Non_beta_simulations.R"); source("_tmp/Heterogeneous_effect_sim.R")'
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
--------------------

--- sbatch Command ---
sbatch  -o /scratch/abrmoe030/projects/mimosa2/_tmp/log/sbatch/sim/run_2026-08-20_09-54-59_out.txt /scratch/abrmoe030/projects/mimosa2/scripts/slurm/sim.sh
-----------------------

--- Standard Output ---
Submitted batch job 1248643
-----------------------

