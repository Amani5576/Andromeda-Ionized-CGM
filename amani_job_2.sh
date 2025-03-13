#!/bin/bash
#SBATCH --job-name="m31_Individual_plots_on_RM_vs_rad_proj"
#SBATCH --output=log_stash/m31-Individual-plots-on-RM-vs-rad-proj-%j.log
#SBATCH --error=log_stash/m31-Individual-plots-on-RM-vs-rad-proj-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=00:30:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for Testing individual plots of Honours |RM| against radial projection"

#Script of interest to be run in current job
singularity exec /idia/software/containers/ASTRO-PY3.simg python3 M31_signal_vs_entire_sky.py 

# Wait for job to complete before logging stats
wait

# Capture job statistics after completion
LOG_FILE="log_stash/m31-Individual-plots-on-RM-vs-rad-proj-$SLURM_JOB_ID.log"
echo "----------------------------------------------------" >> $LOG_FILE
echo "SLURM Job Efficiency Report for Job ID: $SLURM_JOB_ID" >> $LOG_FILE
sacct -j $SLURM_JOB_ID --format=JobID,MaxRSS,Elapsed,State >> $LOG_FILE
echo "----------------------------------------------------" >> $LOG_FILE