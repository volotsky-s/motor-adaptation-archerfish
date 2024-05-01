
# Motor Adaptation in the Archerfish

## Description
This project analyzes the motor adaptation of archerfish to aerial perturbations. 
In our experiments, archerfish were tasked with shooting at an aerial target while we introduced an airflow to perturb their shots. 
The aim was to observe how quickly and effectively the fish could adapt to these disturbances and still hit their targets. 
The experiments were structured into different stages (baseline, perturbation, washout, and in one case, reverse direction perturbation) to analyze the adaptation process.

## Repository Contents
- `motor_adaptation_archerfish.csv`: Main dataset with all trial data: distance (error), trial number, fish number, experiment type, experiment phase, epoch.
- `adaptation_main.R`: Main R script for analysis.
- R function scripts sourced by the main script.
- Output example CSV files: Outputs generated by the R scripts.
- MATLAB scripts (`*.m`): Used for building graphs based on the R-generated CSV files.

### Prerequisites
You will need R and MATLAB installed on your machine. Additionally, the following R packages are required:
- `jags`
- `rjags`
- `runjags`
