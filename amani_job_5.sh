#!/bin/bash
#SBATCH --job-name="m31_histogram_annuli_analysis-COMBINED-%j"
#SBATCH --output=log_stash/m31_histogram_annuli_analysis-COMBINED-%j.log
#SBATCH --error=log_stash/m31_histogram_annuli_analysis-COMBINED-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50G
#SBATCH --time=0:30:00
#SBATCH --partition=Main
#SBATCH --array=1-3  # Run three job instances

set -e  # Exit immediately if any command fails

# Define the singularity command
SINGULARITY_CMD="singularity exec /idia/software/containers/ASTRO-PY3.10-2024-10-18.sif python3 M31_signal_vs_entire_sky.py --annuli-anal"

# Run different tasks based on Slurm array ID
case $SLURM_ARRAY_TASK_ID in
    1)
        echo "Submitting Slurm Job for plotting Histograms of random patches of the sky, with superimposed histograms"
        $SINGULARITY_CMD --overplot
        ;;
    2)
        echo "Submitting Slurm Job for plotting Histograms just for M31 alone, per individual annulus"
        $SINGULARITY_CMD --m31-annuli-anal
        ;;
    3)
        echo "Submitting Slurm Job for plotting Histograms just for M31 alone, with an Overplot of all annuli"
        $SINGULARITY_CMD --m31-annuli-anal --overplot
        ;;
esac
