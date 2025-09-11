setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/')
library(here)
library(SpatialExperiment)
library(scran)
library(tidyverse)
library(escheR)
library(scater)
library(scattermore)

spe <- readRDS(here('processed-data/xenium/raw_combined_spe.RDS'))

plot_coldata_on_tissue <- function(x, column_name){
    plist <- list()
    brnums <- unique(x$sample_type)
    for (i in 1:length(brnums)){
        x_sub <- x[, x$sample_type == brnums[i]]
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

pdf(here("plots", "xenium", "01_QC_metrics.pdf"))
plotColData(spe, y="total_counts", x="sample_type")+geom_scattermore()
plots

plotColData(spe, y="control_probe_counts", x="sample_type")+geom_scattermore()
plots1

plotColData(spe, y="unassigned_codeword_counts", x="sample_type")+geom_scattermore()
plots2

plotColData(spe, y="cell_area", x="sample_type")+geom_scattermore()
plots3

plotColData(spe, y="nucleus_area", x="sample_type")+geom_scattermore()
plots4

plotColData(spe, y="transcript_counts", x="sample_type")+ geom_scattermore()
plots5
dev.off()
