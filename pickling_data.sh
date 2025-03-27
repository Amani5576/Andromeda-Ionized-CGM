#!/bin/bash
#SBATCH --job-name="Overwriting-pickled-data-%j"
#SBATCH --output=log_stash/Overwriting-pickled-data-%j.log
#SBATCH --error=log_stash/Overwriting-pickled-data-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=30G
#SBATCH --time=00:45:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for M31 vs 10K Random Patches Storage and overwriting of Pickled data."

# Script of interest to be run in current job
singularity exec /idia/software/containers/ASTRO-PY3.10-2024-10-18.sif python3 M31_signal_vs_entire_sky.py --pickling