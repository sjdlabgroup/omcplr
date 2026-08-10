#!/usr/bin/env bash

# ==============================================================================
# Pipeline Name: Evolutionary Dynamics & Clustered Genomic Alterations Simulation
# Description: Executes a complete simulation run of clonally amplified genomic 
#              alterations under distinct constraints, infers complex clustered 
#              events with configurable noise and dropout factors, and generates 
#              summarized annotation reports for downstream tools (e.g., omcplr).
#
# Environment: Python 3.9.6+ (Compatible with other Python 3.x versions)
# Dependencies: pandas, numpy
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -euo pipefail

# ==============================================================================
# Configuration Variables & File Paths
# ==============================================================================
SIM_OUTPUT_DIR="tmp/sim_run_30"
INTERMEDIATE_REPORT="output_report.tsv"
FINAL_SUMMARY_REPORT="summary_report.tsv"

# Clustered events inference parameters
MIN_CELLS=1
NOISE_FACTOR=0.2
DROPOUT_FACTOR=0.05

# ==============================================================================
# Pipeline Execution
# ==============================================================================

echo "======================================================================"
echo " [Step 1/3] Running Clone Dynamics Simulation"
echo " Inputs: chromosomes.bed, parameters.txt"
echo "======================================================================"
python3 simulation_clone_dynamics.py


echo -e "\n======================================================================"
echo " [Step 2/3] Identifying Clustered & Overlapping Events"
echo " Target Folder: $SIM_OUTPUT_DIR"
echo " Noise Factor: $NOISE_FACTOR | Dropout Factor: $DROPOUT_FACTOR"
echo "======================================================================"
python3 simulation_clustered_events_inference.py \
    "$SIM_OUTPUT_DIR" \
    "$INTERMEDIATE_REPORT" \
    "$MIN_CELLS" \
    "$NOISE_FACTOR" \
    "$DROPOUT_FACTOR"


echo -e "\n======================================================================"
echo " [Step 3/3] Generating Summarized Annotations"
echo " Input Report: $INTERMEDIATE_REPORT"
echo "======================================================================"
python3 simulation_summary_wrap.py \
    "$INTERMEDIATE_REPORT" \
    "$FINAL_SUMMARY_REPORT"


echo -e "\n======================================================================"
echo " Pipeline Execution Completed Successfully!"
echo " Final summary report successfully written to: $FINAL_SUMMARY_REPORT"
echo "======================================================================"