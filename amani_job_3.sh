#!/bin/bash
#SBATCH --job-name="m31_histogram_annuli_analysis-%j"
#SBATCH --output=log_stash/m31_histogram_annuli_analysis-%j.log
#SBATCH --error=log_stash/m31_histogram_annuli_analysis-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=58G
#SBATCH --time=01:30:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for plotting Histograms after presenting stacked (non-mean) verison at 12th March Wednesday meeting"

#Script of interest to be run in current job
singularity exec /idia/software/containers/ASTRO-PY3.10-2024-10-18.sif python3 M31_signal_vs_entire_sky.py --annuli-anal
