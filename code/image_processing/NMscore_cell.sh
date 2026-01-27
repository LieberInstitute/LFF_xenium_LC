#!/bin/bash
#SBATCH --mem=100G
#SBATCH --job-name=NMscore_cell
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/NMscore_cell_%a.txt
#SBATCH -e logs/NMscore_cell_%a.txt
#SBATCH --array=1-35%5

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load matlab
module list

# --------- CONFIG ---------

brnum=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $1}' inputs.txt | awk '{print $1}')
echo "$brnum"

matlab -nodesktop -nosplash -r "NMscore_cell('$brnum')"

echo "**** Job ends ****"
date