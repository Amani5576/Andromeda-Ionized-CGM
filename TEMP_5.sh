#!/bin/bash

# TEMP_5.sh - Runs the 4 TEMP_5-related scripts in sequence

# Make sure each script is executable
chmod +x temp_5_dispersion.sh
chmod +x temp_5_not_smoothed.sh
chmod +x temp_5.sh
chmod +x temp_5_smoothed.sh

# Run each script
./temp_5_dispersion.sh
./temp_5_not_smoothed.sh
./temp_5.sh
./temp_5_smoothed.sh
