#!/bin/bash
#SBATCH --mem=285G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -o 01_xeniumranger_NM_DAPI_%a.txt
#SBATCH -e 01_xeniumranger_NM_DAPI_%a.txt
#SBATCH --array=5-8%5
#SBATCH --constraint="intel"

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module --ignore_cache load xeniumranger/3.1.1

BASE="/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC"

# Get brnum from Nth line (tab-delimited). Stops after first print.
brnum=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $1}' ../inputs.txt | awk '{print $1}')
echo "brnum = $brnum"

bundle=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'BEGIN{FS="\t"} NR==n {print $6}' ../inputs.txt | awk '{print $1}')
echo "bundle = $bundle"


MASK=${BASE}/processed-data/xenium_imageProcessing/${brnum}/combined_nucmask.npy
BUNDLE=${BASE}/raw-data/xenium/xenium-instrument/${bundle}
OUTPATH=xeniumranger_NM_DAPI_${brnum}

xeniumranger import-segmentation --xenium-bundle=${BUNDLE} --id=${OUTPATH} --nuclei=${MASK} --localcores=48 --localmem=250 --disable-ui=true --jobmode=local

# move output where it should've been able to go in the first place
mv ${OUTPATH} ${BASE}/processed-data/xenium_imageProcessing/${brnum}/
echo "**** Job ends ****"
date
