#!/bin/bash
#SBATCH --job-name="m31_Individual_plots_on_RM_vs_rad_proj"
#SBATCH --output=log_stash/m31-Individual-plots-on-RM-vs-rad-proj-%j.log
#SBATCH --error=log_stash/m31-Individual-plots-on-RM-vs-rad-proj-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=35G
#SBATCH --time=00:40:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for Testing individual plots of Honours |RM| against radial projection"

#Script of interest to be run in current job
singularity exec /idia/software/containers/ASTRO-PY3.simg python3 M31_signal_vs_entire_sky.py