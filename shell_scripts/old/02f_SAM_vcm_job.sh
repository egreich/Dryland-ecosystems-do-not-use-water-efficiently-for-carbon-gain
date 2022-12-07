#!/bin/bash
#SBATCH --job-name=SAM_vcm
##SBATCH --workdir
#SBATCH --cpus-per-task=8
#SBATCH --time=30:00:00
#SBATCH --mem=150000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu

chmod +x ./scripts/02f_script_SAM_vcm.R
srun ./scripts/02f_script_SAM_vcm.R
