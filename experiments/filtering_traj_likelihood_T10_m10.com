#!/bin/bash
#SBATCH --job-name=trajectory_T10_m10
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}/.."

Rscript experiments/filtering_traj_likelihood_T10_m10.R
