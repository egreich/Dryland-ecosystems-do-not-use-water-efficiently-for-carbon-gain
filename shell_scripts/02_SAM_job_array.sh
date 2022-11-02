#!/bin/bash
#SBATCH --job-name=SAM_$a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/log/AllSites_%A_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=150000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=1-21

### %A is monsoon job number %a is interior array index

chain=$(sed -n "$SLURM_ARRAY_TASK_ID"p chainEND)
site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedEND)

chmod +x ./scripts/02_script_SAM_HPC.R # for permissions
srun ./shell_scripts/run_SAM_job.sh $chain $site $seed
