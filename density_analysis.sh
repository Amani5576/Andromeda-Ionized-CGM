#!/bin/bash
for bw in $(seq 0.01 0.005 0.1); do python3 M31_signal_density.py --bw $bw --save-plot --cumulative; done
