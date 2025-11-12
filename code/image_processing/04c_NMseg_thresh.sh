#!/bin/bash
#SBATCH --mem=60G
#SBATCH --job-name=NMsegThresh
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/NMseg_thresh_%a.txt
#SBATCH -e logs/NMseg_thresh_%a.txt
#SBATCH --array=2,3,5,11-13,15,17,28,30%5

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
K=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $4}' inputs.txt | awk '{print $1}')
echo "$K"

matlab -nodesktop -nosplash -r "NMseg_thresh('$brnum',$K)"

echo "**** Job ends ****"
date