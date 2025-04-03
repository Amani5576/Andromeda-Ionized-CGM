#!/bin/bash
#SBATCH --job-name="Pickle-and-Original-Plot-2order-BG-Subtraction%j"
#SBATCH --output=log_stash/Pickle-and-Original-Plot-2order-BG-Subtraction-%j.log
#SBATCH --error=log_stash/Pickle-and-Original-Plot-2order-BG-Subtraction-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=07:00:00
#SBATCH --partition=Main

echo "Submitting Slurm Job: BG Subtraction with Spline, Pickling, and Plotting."

# Run the Python script for BG subtraction, saving plots, and pickling data
singularity exec /idia/software/containers/ASTRO-PY3.10-2024-10-18.sif python3 M31_signal_vs_entire_sky.py --original-plot --save-plot --pickling

echo "Job completed: Pickling and Original Plot generation finished."
