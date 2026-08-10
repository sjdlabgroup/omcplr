# Simulation Framework for Genomic Rearrangement Inference

## Description
# This simulation framework is designed to model the evolutionary dynamics of clonally amplified genomic alterations. It consists of three modular Python scripts and a master shell script that coordinates the pipeline.

#The framework facilitates the study of:
# Simulates stochastic cell birth, death, DNA breaks, and the accumulation of rearrangements across subclonal lineages.
# Identifies genomic regions with overlapping rearrangements, extracting individual and combinatorial count data suitable for downstream tools such as `omcplr`.
# Allows for the fine-tuning of genomic instability, catastrophic event complexity, and the introduction of technical noise and dropout rates to assess model robustness.

## Prerequisites
# Python 3.x
# Standard scientific libraries (e.g., `numpy`, `pandas`)

## Usage
# The master script executes a complete simulation run, from initial population dynamics to the generation of summarized annotation reports.

# ```bash
# sh run_simulation_complex_event_generation.sh
