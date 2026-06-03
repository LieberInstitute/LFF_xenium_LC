setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/')
library(here)
library(SpatialExperiment)
library(SingleCellExperiment)
library(readxl)
library(tidyverse)


########################################################
# Read in the raw Xenium data into SPE objects and save.
#######################################################
brnum = Sys.getenv("BRNUM")
sample_path = here('processed-data','xenium_imageProcessing_new/xeniumRanger/',paste0('xeniumranger_NM_DAPI_', brnum), 'outs')

    counts_path <- here(sample_path, "cell_feature_matrix.h5")
    cell_info_path <- here(sample_path, "cells.csv.gz")

    sce <- DropletUtils::read10xCounts(counts_path)
    counts(sce) <- methods::as(DelayedArray::realize(counts(sce)), "dgCMatrix") # Convert to delayed array

    cell_info <- vroom::vroom(cell_info_path)

    colData(sce) <- cbind(colData(sce), cell_info)
    spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x_centroid", "y_centroid"))
    rownames(spe) <- rowData(spe)$Symbol # change rownames to gene symbol

   
saveRDS(spe, here("processed-data/01_spe/NMDAPInewseg_rawSPE", paste0(brnum,".RDS")))
