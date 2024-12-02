#!/bin/bash
#SBATCH --job-name=RUN_FAETH
#SBATCH --time=20:00:00
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --mem=40000
#SBATCH --ntasks=4
#SBATCH --account=HG
#SBATCH --qos=std

module load gcc
module load R/4.0.2
module load SHARED/calc_grm/main

# setup R library path
R_LIBS_USER="home/HG/perez02/R/x86_64-pc-linux-gnu-library/4.0"

#store the base path (main folder)
BASE_PATH=$(dirname "$(pwd)")

# construct the path to GNSW_FAETH_Pipeline relative to the current directory
SCRIPT_PATH="$BASE_PATH/1_Scripts/GNSW_FAETH_VCE.R"

# run StartScript.R - stert the simulation pipeline
Rscript "$SCRIPT_PATH" "$BASE_PATH" > VCE.log