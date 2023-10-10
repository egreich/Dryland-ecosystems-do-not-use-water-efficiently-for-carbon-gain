#!/bin/bash
#SBATCH --job-name=relinf_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/SAM/log/reinf_%A_%a.log
##SBATCH --cpus-per-task=3
#SBATCH --time=40:00:00
#SBATCH --mem=80000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
#SBATCH --array=1-48

### %A is monsoon job number %a is interior array index

module load R/4.1.2 # load a specific R version

chmod +x ./scripts/03_script_relinf.R # for permissions
chmod +x ./shell_scripts/run_relinf_job.sh # for permissions


site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteinfEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedinfEND)
modelv=$(sed -n "$SLURM_ARRAY_TASK_ID"p modelvinfEND)
voi=$(sed -n "$SLURM_ARRAY_TASK_ID"p voiEND)

# Run the analysis
srun ./shell_scripts/run_relinf_job.sh $site $seed $modelv $voi
