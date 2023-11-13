#!/bin/bash
#SBATCH --job-name=postSAM_%a
##SBATCH --workdir
#SBATCH --output=/scratch/egr65/SAM/log/postSAM_%A_%a.log
##SBATCH --cpus-per-task=3
#SBATCH --time=2:00:00
#SBATCH --mem=50000
#SBATCH --mail-type=all
#SBATCH --mail-user=egr65@nau.edu
##SBATCH --array=10,11,12,13,17,18,19,20,21,23,9,33-36
##SBATCH --array=already ran 1-8,14-16,22,24-27,29-32,38-40
#SBATCH --array=1-40
## not 10, 11, 12, 13, 17,18,19,20,21,23,9, 33-36
## 28 problem with sig.WUE site4 model 8
## 37 Compilation error on line 263. Possible directed cycle involving var.pred. site 5 model 9

### %A is monsoon job number %a is interior array index

module load R/4.1.2 # load a specific R version

chmod +x ./scripts/02_script_SAM.R # for permissions
chmod +x ./shell_scripts/run_SAM_job.sh # for permissions

site=$(sed -n "$SLURM_ARRAY_TASK_ID"p siteEND)
seed=$(sed -n "$SLURM_ARRAY_TASK_ID"p seedEND)
modelv=$(sed -n "$SLURM_ARRAY_TASK_ID"p modelvEND)

# Run the analysis
srun ./shell_scripts/run_SAM_job.sh $site $seed $modelv
