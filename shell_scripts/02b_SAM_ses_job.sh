#!/bin/bash
#SBATCH --job-name=SAM_ses
##SBATCH --workdir
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=150000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu

chmod +x ./scripts/02b_script_SAM_ses.R
srun ./scripts/02b_script_SAM_ses.R
