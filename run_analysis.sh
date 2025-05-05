#!/bin/bash

# Bash script to run a sequence of analysis scripts

python3 M31_signal_vs_entire_sky.py --pickling

python3 M31_signal_vs_entire_sky.py --rm-vs-PA

python3 M31_signal_vs_entire_sky.py --rm-vs-gal-lat

python3 RM_vs_M31_radius.py --save-plot

python3 M31_signal_density.py --save-plot

python3 M31_signal_density.py --cumulative --save-plot

echo "All scripts run by 'run_analysis' executed successfully."
