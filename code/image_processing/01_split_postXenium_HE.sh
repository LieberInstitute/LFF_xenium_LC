#!/bin/bash
#SBATCH --mem=40G
#SBATCH --job-name=split_postXeniumHE
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/split_postXeniumHE_%a.txt
#SBATCH -e logs/split_postXeniumHE_%a.txt
#SBATCH --array=1-29%5

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
SLIDES_DIR="/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/raw-data/xenium/post-xenium_images"
slides=( $(ls -1 "$SLIDES_DIR") )
idx=$((SLURM_ARRAY_TASK_ID - 1))
sample="${slides[$idx]}"    

echo "Processing:  ${sample}"
python 01_split_postXenium_HE.py "${sample}"

echo "**** Job ends ****"
date