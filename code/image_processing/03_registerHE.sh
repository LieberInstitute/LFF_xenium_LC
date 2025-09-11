#!/bin/bash
#SBATCH --mem=120G
#SBATCH --job-name=registerHE
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/registerHE_%a.txt
#SBATCH -e logs/registerHE_%a.txt
#SBATCH --array=5

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load visium_hd/1.0
module list

# --------- CONFIG ---------

brnum=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $1}' inputs.txt | awk '{print $1}')
echo "$brnum"
slide=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $2}' inputs.txt | awk '{print $1}')
echo "$slide"
sample=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $3}' inputs.txt | awk '{print $1}')
echo "$sample"

python 03_registerHE.py $brnum $slide $sample

echo "**** Job ends ****"
date