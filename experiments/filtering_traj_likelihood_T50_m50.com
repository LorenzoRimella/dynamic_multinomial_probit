#!/bin/bash
#SBATCH --job-name=trajectory_T50_m50
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/.."

Rscript experiments/filtering_traj_likelihood_T50_m50.R
