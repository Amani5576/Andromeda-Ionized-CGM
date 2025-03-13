#!/bin/bash
#SBATCH --job-name="m31_Individual_plots_on_RM_vs_rad_proj"
#SBATCH --output=log_stash/m31-Individual-plots-on-RM-vs-rad-proj-%j.log
#SBATCH --error=log_stash/m31-Individual-plots-on-RM-vs-rad-proj-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:18:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for plotting Histograms after presenting stacked (non-mean) verison at 12th March Wednesday meeting"

#Script of interest to be run in current job
singularity exec /idia/software/containers/ASTRO-PY3.simg python3 M31_signal_vs_entire_sky.py --annuli-anal