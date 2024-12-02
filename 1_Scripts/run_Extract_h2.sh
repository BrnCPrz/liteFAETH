#!/bin/bash
#SBATCH --job-name=RUN_Exth2
#SBATCH --time=10:00:00
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --mem=10000
#SBATCH --ntasks=2
#SBATCH --account=HG
#SBATCH --qos=std

module load gcc
module load R/4.0.2
module load SHARED/calc_grm/main

# setup R library path
R_LIBS_USER="home/HG/perez02/R/x86_64-pc-linux-gnu-library/4.0"

#store the path to working directory
BASE_PATH="$PWD"

# construct the path to GNSW_FAETH_Extract_h2.R relative to the current directory
SCRIPT_PATH="$BASE_PATH/GNSW_FAETH_Extract_h2.R"

# run StartScript.R - stert the simulation pipeline
Rscript "$SCRIPT_PATH" "$BASE_PATH" > Extract_h2.log

