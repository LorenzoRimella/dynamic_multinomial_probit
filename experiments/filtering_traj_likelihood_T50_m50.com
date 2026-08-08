#!/bin/bash
#SBATCH --job-name=trajectory_T50_m50
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=g100_usr_pmem      # Partition (queue) to submit to
#SBATCH --qos=g100_qos_lprod                  # Quality of Service
#SBATCH --mem=64G                          # Memory per node (adjust as needed)
#SBATCH --account=uBG25_EcoRim       # Project account number

# Load necessary modules
module load profile/bioinf
module load autoload r

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/.."

Rscript experiments/filtering_traj_likelihood_T50_m50.R
