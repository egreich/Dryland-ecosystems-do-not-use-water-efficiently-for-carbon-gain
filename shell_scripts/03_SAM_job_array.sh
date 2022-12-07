#!/bin/bash
#SBATCH --job-name=SAM_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/SAM/log/AllSites_%A_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --time=90:00:00
#SBATCH --mem=150000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=1-24

### %A is monsoon job number %a is interior array index

chmod +x ./scripts/03_script_SAM_HPC.R # for permissions

chain=$(sed -n "$SLURM_ARRAY_TASK_ID"p chainEND)
site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedEND)

# Run the analysis
srun ./shell_scripts/run_SAM_job.sh $chain $site $seed
