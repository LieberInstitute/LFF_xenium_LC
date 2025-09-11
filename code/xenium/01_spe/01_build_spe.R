setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/')
library(here)
library(SpatialExperiment)
library(SingleCellExperiment)
library(readxl)
library(tidyverse)


########################################################
# Read in the raw Xenium data into SPE objects and save.
#######################################################

sample_path = 'raw-data/xenium/xenium-instrument/output-XETG00558__0068654__Br6297__20250501__172909'

    counts_path <- here(sample_path, "cell_feature_matrix.h5")
    cell_info_path <- here(sample_path, "cells.csv.gz")

    sce <- DropletUtils::read10xCounts(counts_path)
    counts(sce) <- methods::as(DelayedArray::realize(counts(sce)), "dgCMatrix") # Convert to delayed array

    cell_info <- vroom::vroom(cell_info_path)

    colData(sce) <- cbind(colData(sce), cell_info)
    spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x_centroid", "y_centroid"))
    rownames(spe) <- rowData(spe)$Symbol # change rownames to gene symbol

   
saveRDS(spe, here("processed-data/xenium/DAPI", "raw_spe.RDS"))