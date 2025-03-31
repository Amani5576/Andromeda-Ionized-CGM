#!/bin/bash
#SBATCH --job-name="Overwriting-pickled-data-BG-Subtraction-Spline-%j"
#SBATCH --output=log_stash/Overwriting-pickled-data-BG-Subtraction-Spline-%j.log
#SBATCH --error=log_stash/Overwriting-pickled-data-BG-Subtraction-Spline-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=38G
#SBATCH --time=01:30:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job: Proper BG subtraction with spline and pickling."

# Run the Python script for BG subtraction with pickling enabled
singularity exec /idia/software/containers/ASTRO-PY3.10-2024-10-18.sif python3 M31_signal_vs_entire_sky.py --pickling