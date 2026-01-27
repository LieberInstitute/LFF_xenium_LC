#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=32G
#SBATCH --job-name=parquet2csv
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o logs/parquet2csv.log
#SBATCH -e logs/parquet2csv.log

set -euo pipefail

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"

repo_dir='/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC'
inputs="$repo_dir/code/image_processing/inputs.txt"

module load visium_hd/1.0

# Read first column into a bash array (split on whitespace/newlines)
mapfile -t SAMPLES < <(awk -F'\t' '{print $1}' "$inputs")

for SAMPLE in "${SAMPLES[@]}"; do
  spatial_dir="$repo_dir/processed-data/xenium_imageProcessing/${SAMPLE}/xeniumranger_NM_DAPI_${SAMPLE}/outs"

  echo "---- ${SAMPLE} ----"
  echo "outs: $spatial_dir"

  if [[ ! -d "$spatial_dir" ]]; then
    echo "WARNING: outs dir not found, skipping: $spatial_dir"
    continue
  fi

  # Convert cells.parquet -> cells.csv if needed
  if [[ -f "$spatial_dir/cells.parquet" && ! -f "$spatial_dir/cells.csv" ]]; then
    echo "Converting cells.parquet -> cells.csv"
    parquet-tools csv "$spatial_dir/cells.parquet" > "$spatial_dir/cells.csv"
  fi

  # Convert cell_boundaries.parquet -> cell_boundaries.csv if needed
  if [[ -f "$spatial_dir/cell_boundaries.parquet" && ! -f "$spatial_dir/cell_boundaries.csv" ]]; then
    echo "Converting cell_boundaries.parquet -> cell_boundaries.csv"
    parquet-tools csv "$spatial_dir/cell_boundaries.parquet" > "$spatial_dir/cell_boundaries.csv"
  fi
done

echo "**** Job ends ****"
date
