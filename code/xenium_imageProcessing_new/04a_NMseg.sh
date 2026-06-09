#!/bin/bash
#SBATCH --job-name=NMseg
#SBATCH --mem=80G
#SBATCH -o logs/NMseg_%a.txt
#SBATCH -e logs/NMseg_%a.txt
#SBATCH --array=1-33%10

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load MATLAB
module load matlab

brnum=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $1}' inputs.txt | awk '{print $1}')
echo "$brnum"

matlab -nodesktop -nosplash -r "NMseg('$brnum')"

echo "**** Job ends ****"
date