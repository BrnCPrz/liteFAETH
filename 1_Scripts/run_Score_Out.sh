#!/bin/bash
#SBATCH --job-name=RUN_ScrOut
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

# construct the path to relative to the current directory
SCRIPT_PATH="$BASE_PATH/1_Scripts/GNSW_FAETH_Score_Out.R"

# run Script.R 
Rscript "$SCRIPT_PATH" "$BASE_PATH" > Score_Out.log