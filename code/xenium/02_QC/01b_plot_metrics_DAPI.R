setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/')
library(here)
library(SpatialExperiment)
library(scran)
library(tidyverse)
library(escheR)
library(scater)
library(scattermore)

spe <- readRDS(here('processed-data/xenium/DAPI_rawSPE/raw_combined_spe.RDS'))

plot_coldata_on_tissue <- function(x, column_name){
    plist <- list()
    brnums <- unique(x$brnum)
    for (i in 1:length(brnums)){
        x_sub <- x[, x$brnum == brnums[i]]
         p <- (make_escheR(x_sub) %>%
            add_fill(column_name)+
            ggtitle(brnums[i])+
            geom_scattermore()
			)
			 plist[[i]] = p
    }
    return(plist)
}

plots <- plot_coldata_on_tissue(spe, "total_counts")
plots1 <- plot_coldata_on_tissue(spe, "control_probe_counts")
plots2 <- plot_coldata_on_tissue(spe, "unassigned_codeword_counts")
plots3 <- plot_coldata_on_tissue(spe, "cell_area")
plots4 <- plot_coldata_on_tissue(spe, "nucleus_area")
plots5 <- plot_coldata_on_tissue(spe, "transcript_counts")

pdf(here("plots", "xenium1", "02_QC", "total_counts_DAPI.pdf"))
plots
dev.off()

pdf(here("plots", "xenium1", "02_QC", "control_probe_counts_DAPI.pdf"))
plots1
dev.off()

pdf(here("plots", "xenium1", "02_QC", "unassigned_codeword_counts_DAPI.pdf"))
plots2
dev.off()

pdf(here("plots", "xenium1", "02_QC", "cell_area_DAPI.pdf"))
plots3
dev.off()

pdf(here("plots", "xenium1", "02_QC", "nucleus_area_DAPI.pdf"))
plots4
dev.off()

pdf(here("plots", "xenium1", "02_QC", "transcript_counts_DAPI.pdf"))
plots5
dev.off()
