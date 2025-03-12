#!/bin/bash
#SBATCH --job-name="m31_Individual_plots_on_RM_vs_rad_proj"
#SBATCH --output=m31-Individual-plots-on-RM-vs-rad-proj-%j.log
#SBATCH --error=m31-Individual-plots-on-RM-vs-rad-proj-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=00:50:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for Testing individual plots of Honours |RM| against radial projection"

# Running my script
singularity exec /idia/software/containers/ASTRO-PY3.simg python3 M31_signal_vs_entire_sky.py 

# Checking if the script executed successfully
if [ $? -eq 0 ]; then
    echo "run successfully"
else
    echo "There was an issue. Run lerr in terminal for more info"
