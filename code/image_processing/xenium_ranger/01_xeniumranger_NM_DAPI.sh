#!/bin/bash
#SBATCH --mem=285G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH -o /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/code/xenium_imageProcessing/xenium_ranger/01_xeniumranger_NM_DAPI.txt
#SBATCH -e /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/code/xenium_imageProcessing/xenium_ranger/01_xeniumranger_NM_DAPI.txt
#SBATCH --array=1
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

# get current sample to resegment
#FILE=/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/raw-data/xenium/output-XETG00558__0068654__Br6297__20250501__172909
FILE=/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/_xenium_LIBD/250505_run-2/output-XETG00558__0068654__Br6297__20250501__172909
MASK=/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/Br6297/combined_nucmask.npy
OUTPATH=xeniumranger_NM_DAPI

xeniumranger import-segmentation --xenium-bundle=${FILE} --id=${OUTPATH} --nuclei=${MASK} --localcores=48 --localmem=250 --disable-ui=true --jobmode=local

# move output where it should've been able to go in the first place
mv ${OUTPATH} /dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/Br6297/
echo "**** Job ends ****"
date
