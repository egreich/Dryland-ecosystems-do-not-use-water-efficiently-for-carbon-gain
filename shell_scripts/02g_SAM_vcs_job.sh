#!/bin/bash
#SBATCH --job-name=SAM_vcs
##SBATCH --workdir
#SBATCH --cpus-per-task=8
#SBATCH --time=20:00:00
#SBATCH --mem=150000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu

chmod +x ./scripts/02g_script_SAM_vcs.R
srun ./scripts/02g_script_SAM_vcs.R
