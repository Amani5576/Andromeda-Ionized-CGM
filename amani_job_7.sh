#!/bin/bash
#SBATCH --job-name="m31-vs-10K-patches-%j"
#SBATCH --output=log_stash/m31-vs-10K-patches-%j.log
#SBATCH --error=log_stash/m31-vs-10K-patches-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for M31 vs 10K Random Patches Analysis. Counting how many RM landed in each patch"
echo "Plotting Total RM values per patch."

# Script of interest to be run in current job
singularity exec /idia/software/containers/ASTRO-PY3.10-2024-10-18.sif python3 M31_signal_vs_entire_sky.py --rm-per-patch-hist
