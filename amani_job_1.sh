!/bin/bash
#SBATCH --job-name="AndromedaAnnulusCalc"
#SBATCH --output=Andromeda-annulus-calc-%j.log
#SBATCH --error=Andromeda-annulus-calc-%j-error.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=00:50:00
#SBATCH --partition=Main

# Load required modules
echo "Submitting Slurm Job for Andromeda Annulus Calculations"

# Running my script
singularity exec /idia/software/containers/ASTRO-PY3.simg python3 M31_signal_vs_entire_sky.py 

echo "Done Running"