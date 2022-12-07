#!/bin/bash
#SBATCH --job-name=onesite_SAM_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/SAM/log/OneSite_%A_%a.log
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=100000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=1-3

### %A is monsoon job number %a is interior array index

chmod +x ./shell_scripts/run_SAM_job.sh # for permissions
chmod +x ./scripts/03_script_SAM_HPC.R # for permissions

chain=$(sed -n "$SLURM_ARRAY_TASK_ID"p chainoneEND)
site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteoneEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedoneEND)

# Run the analysis
srun ./shell_scripts/run_SAM_job.sh $chain $site $seed
