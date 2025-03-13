#!/bin/bash
#SBATCH --job-name="m31_Individual_plots_on_RM_vs_rad_proj"
#SBATCH --output=log_stash/m31-Individual-plots-on-RM-vs-rad-proj-%j-$(date +%Y%m%d-%H%M%S).log
#SBATCH --error=log_stash/m31-Individual-plots-on-RM-vs-rad-proj-%j-$(date +%Y%m%d-%H%M%S)-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=00:30:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for Testing individual plots of Honours |RM| against radial projection"

# Running main.py first to define variables
singularity exec /idia/software/containers/ASTRO-PY3.simg python3 main.py

#Script of interest to be run in current job
singularity exec /idia/software/containers/ASTRO-PY3.simg python3 M31_signal_vs_entire_sky.py 


# Capture job statistics after completion
echo "----------------------------------------------------" >> m31-Individual-plots-on-RM-vs-rad-proj-$SLURM_JOB_ID.log
echo "SLURM Job Efficiency Report for Job ID: $SLURM_JOB_ID" >> m31-Individual-plots-on-RM-vs-rad-proj-$SLURM_JOB_ID.log
seff $SLURM_JOB_ID >> m31-Individual-plots-on-RM-vs-rad-proj-$SLURM_JOB_ID.log
echo "----------------------------------------------------" >> m31-Individual-plots-on-RM-vs-rad-proj-$SLURM_JOB_ID.log
