#!/bin/bash
#SBATCH --job-name="m31_histogram_annuli_analysis-%j"
#SBATCH --output=log_stash/m31_histogram_annuli_analysis-%j.log
#SBATCH --error=log_stash/m31_histogram_annuli_analysis-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=00:30:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for plotting Histograms after presenting stacked (non-mean) verison at 12th March Wednesday meeting"

#Script of interest to be run in current job
singularity exec /idia/software/containers/casa-stable.img casa --log2term --nologger
