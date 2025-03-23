#!/bin/bash
#SBATCH --job-name="m31-vs-10K-random-patches-KS-Test-%j"
#SBATCH --output=log_stash/m31-vs-10K-random-patches-KS-Test-%j.log
#SBATCH --error=log_stash/m31-vs-10K-random-patches-KS-Test-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=30G
#SBATCH --time=1:00:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for KS Testing of M31 vs 10 thousand Virial Radii analogous to M31"
echo "Note: BG correction in each patch has been acheived for all raw RM, from patch to patch"

#Script of interest to be run in current job
singularity exec /idia/software/containers/ASTRO-PY3.10-2024-10-18.sif python3 M31_signal_vs_entire_sky.py --m31-ks-test