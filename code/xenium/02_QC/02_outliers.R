setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/')
library(tidyverse)
library(escheR)
library(here)
library(SpatialExperiment)
library(scuttle)
library(scattermore)
library(gridExtra)
library(matrixStats)


spe <- readRDS(here('processed-data/xenium/raw_combined_spe.RDS'))
rowData(spe)$Gene <- rownames(rowData(spe))
colnames(spe) <- paste0(spe$Barcode, "_", spe$sample_id)
#gene_list <- readxl::read_excel(here("raw-data", "experiment_info", 
#    "Xenium_SCZ_ProbeSelection5_SHK_v4_15_w_ensemblgene_list.xlsx"))
#rowData(spe) <- merge(gene_list, rowData(spe), by="Gene",  all.y=TRUE)


is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
is_mt <- stringr::str_detect(rownames(spe), "MT-")

sample_types <- unique(spe$sample_type)

# Remove empty cells
empty_cells <- colnames(spe)[colSums(counts(spe)) == 0]
spe <- spe[, colSums(counts(spe)) > 0]

all_outlier_ids <- c(empty_cells)
all_outlier_ids
for(i in 1:length(sample_types)){
    br_use <- sample_types[[i]]
    print(sprintf("------------%s------------", br_use))
    spe_sub <- spe[, colData(spe)$sample_type == br_use]
    # Compute within-sample QC metrics based on negative controls and mitochondrial genes
    spe_sub <- scuttle::addPerCellQCMetrics(spe_sub, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        mito=is_mt))
                                


    # Identify the outliers based on these QC metrics
    # For the negative controls, the median is always 0, so we will use a fixed threshold
    neg_probe_quantile <- quantile(spe_sub$subsets_negProbe_percent, 0.99, na.rm=TRUE)
    print(neg_probe_quantile)
    spe_sub$neg_probe_out <- spe_sub$subsets_negProbe_percent >= neg_probe_quantile

    neg_codeword_quantile <- quantile(spe_sub$subsets_negCodeword_percent, 0.99, na.rm=TRUE)
    print(neg_codeword_quantile)
    spe_sub$neg_codeword_out <- spe_sub$subsets_negCodeword_percent >= neg_codeword_quantile

    unassigned_quantile <- quantile(spe_sub$subsets_unassigned_percent, 0.99, na.rm=TRUE) # why are there NAs?
    print(unassigned_quantile)
    spe_sub$unassigned_out <- spe_sub$subsets_unassigned_percent >= unassigned_quantile

    spe_sub$mito_out <- isOutlier(spe_sub$subsets_mito_percent, type="higher")
    spe_sub$detected_out <- isOutlier(spe_sub$detected, type="lower")
    spe_sub$total_counts_out <- isOutlier(spe_sub$total_counts)


    print("any NAs in unassigned_percent?")
    print(any(is.na(spe_sub$subsets_unassigned_percent)))
    print("any NAs in neg_codeword?")
    print(any(is.na(spe_sub$subsets_negCodeword_percent)))


    print("Proportion of negative probe outliers:")
    print(sum(spe_sub$neg_probe_out, na.rm=TRUE)/dim(spe_sub)[[2]])
    print("Proportion of unassigned outliers:")
    print(sum(spe_sub$unassigned_out, na.rm=TRUE)/dim(spe_sub)[[2]])
    print("Proportion of neagtive codeword outliers:")
    print(sum(spe_sub$neg_codeword_out, na.rm=TRUE)/dim(spe_sub)[[2]])
    print("Proportion of mito outliers:")
    print(sum(spe_sub$mito_out, na.rm=TRUE)/dim(spe_sub)[[2]])
    print("Proportion of detected outliers:")
    print(sum(spe_sub$detected_out, na.rm=TRUE)/dim(spe_sub)[[2]])
    print("Proportion of total counts outliers:")
    print(sum(spe_sub$total_counts_out, na.rm=TRUE)/dim(spe_sub)[[2]])


    pdf(here("plots", "xenium", paste0("01_", br_use, ".pdf")), height=15, width=25)
        # Make some plots based on these metrics
    # Plot histogram of subsets_negProbe_percent
    median_negProbe <- median(spe_sub$subsets_negProbe_percent, na.rm=TRUE)
    print(ggplot(colData(spe_sub), aes(x=subsets_negProbe_percent))+
        geom_histogram()+
        geom_vline(aes(xintercept = median_negProbe), color = "red", linetype = "dashed", linewidth = 1) +
        geom_vline(aes(xintercept = neg_probe_quantile), color = "blue", linetype = "dashed", linewidth = 1) +
        ggtitle(br_use))
    
    # Plot histogram of subsets_negCodeword_percent
    median_negCodeword <- median(spe_sub$subsets_negCodeword_percent, na.rm=TRUE)
    print(ggplot(colData(spe_sub), aes(x=subsets_negCodeword_percent))+
        geom_histogram()+
        geom_vline(aes(xintercept = median_negCodeword), color = "red", linetype = "dashed", linewidth = 1) +
        geom_vline(aes(xintercept = neg_codeword_quantile), color = "blue", linetype = "dashed", linewidth = 1) +
        ggtitle(br_use))

    # Plot histogram of subsets_unassigned_percent
    median_unassigned <- median(spe_sub$subsets_negCodeword_percent, na.rm=TRUE)
    print(ggplot(colData(spe_sub), aes(x=subsets_unassigned_percent))+
        geom_histogram()+
        geom_vline(aes(xintercept = median_unassigned), color = "red", linetype = "dashed", linewidth = 1) +
        geom_vline( aes(xintercept = unassigned_quantile), color = "blue", linetype = "dashed", linewidth = 1) +
        ggtitle(br_use))
    
    # Make tissue plots 
    print(make_escheR(spe_sub)%>%
        add_fill("subsets_negProbe_percent")%>%
        add_ground("neg_probe_out")+
        ggtitle(br_use)+
        geom_scattermore())

    print(make_escheR(spe_sub)%>%
        add_fill("subsets_negCodeword_percent") %>%
        add_ground("neg_codeword_out")+
        ggtitle(br_use)+
        geom_scattermore())
    
    print(make_escheR(spe_sub)%>%
        add_fill("subsets_unassigned_percent")%>%
        add_ground("unassigned_out")+
        ggtitle(br_use)+
        geom_scattermore())
    
    print(make_escheR(spe_sub)%>%
        add_fill("subsets_mito_percent")%>%
        add_ground("mito_out")+
        ggtitle(br_use)+
        geom_scattermore())
    
    print(make_escheR(spe_sub)%>%
        add_fill("detected")%>%
        add_ground("detected_out")+
        ggtitle(br_use)+
        geom_scattermore())
    
    print(make_escheR(spe_sub)%>%
        add_fill("total_counts")%>%
        add_ground("total_counts_out")+
        ggtitle(br_use)+
        geom_scattermore())
    
    
    # Make some scatterplots of the QC metrics to see how they are related
    # Inspo for the plots: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009290
    mito_1 <- ggplot(colData(spe_sub), aes(x=detected, y=subsets_mito_percent))+
        geom_point(aes(colour=mito_out))+
        geom_smooth(method="lm")+
        ggtitle(br_use)
    
    mito_2 <- ggplot(colData(spe_sub), aes(x=total_counts, y=subsets_mito_percent))+
        geom_point(aes(colour=mito_out))+
        geom_smooth(method="lm")+
        ggtitle(br_use)

    grid.arrange(mito_1, mito_2, ncol=2)

   negProbe_1 <- ggplot(colData(spe_sub), aes(x=detected, y=subsets_negProbe_percent))+
       geom_point(aes(colour=neg_probe_out))+
       geom_smooth(method="lm")+
       ggtitle(br_use)
   
   negProbe_2 <- ggplot(colData(spe_sub), aes(x=total_counts, y=subsets_negProbe_percent))+
       geom_point(aes(colour=neg_probe_out))+
       geom_smooth(method="lm")+
       ggtitle(br_use)
   
   grid.arrange(negProbe_1, negProbe_2, ncol=2)
  
   negCode_1 <- ggplot(colData(spe_sub), aes(x=detected, y=subsets_negCodeword_percent))+
       geom_point(aes(colour=neg_codeword_out))+
       geom_smooth(method="lm")+
       ggtitle(br_use)
   
   negCode_2 <- ggplot(colData(spe_sub), aes(x=total_counts, y=subsets_negCodeword_percent))+
       geom_point(aes(colour=neg_codeword_out))+
       geom_smooth(method="lm")+
       ggtitle(br_use)
   
   grid.arrange(negCode_1, negCode_2, ncol=2)
  
   unassigned_1 <- ggplot(colData(spe_sub), aes(x=detected, y=subsets_unassigned_percent))+
       geom_point(aes(colour=unassigned_out))+
       geom_smooth(method="lm")+
       ggtitle(br_use)
   
   unassigned_2 <- ggplot(colData(spe_sub), aes(x=total_counts, y=subsets_unassigned_percent))+
       geom_point(aes(colour=unassigned_out))+
       geom_smooth(method="lm")+
       ggtitle(br_use)

    grid.arrange(unassigned_1, unassigned_2, ncol=2)

    rowData(spe_sub)$means <- rowMeans(counts(spe_sub))
    rowData(spe_sub)$vars <- rowVars(counts(spe_sub))

    row_df <- rowData(spe_sub) %>%
        as.data.frame() %>%
        filter(Type=="Gene Expression")


   #mean_probes <- ggplot(row_df, aes(x=Probesets, y=means))+
   #    geom_point()+
   #    ggtitle(br_use)
    
   #var_probes <- ggplot(row_df, aes(x=Probesets, y=vars))+ 
   #    geom_point()+
   #    ggtitle(br_use)
    
    #mean_var <- 
	ggplot(row_df, aes(x=means, y=vars))+
        geom_point()+
        ggtitle(br_use)

    #grid.arrange(mean_var, mean_probes, var_probes, ncol=3)

    # aggregate all outliers by taking the union
    spe_sub$is_outlier <- spe_sub$neg_probe_out | spe_sub$neg_codeword_out | spe_sub$unassigned_out | spe_sub$detected_out | spe_sub$total_counts_out

    print("Proportion of total outliers:")
    print(sum(spe_sub$is_outlier, na.rm=TRUE)/dim(spe_sub)[[2]])
    # tissue plot of all outliers
    print(make_escheR(spe_sub)%>%
        add_fill("is_outlier")+
        ggtitle(br_use)+
        geom_scattermore())
    
    dev.off() 

    outlier_ids <- colnames(spe_sub)[spe_sub$is_outlier] 
    print(length(outlier_ids))   
    all_outlier_ids <- c(all_outlier_ids, outlier_ids)                                               
}

all_outlier_ids <- unique(all_outlier_ids)
print(head(all_outlier_ids))
print(length(all_outlier_ids))
write.csv(all_outlier_ids, here("processed-data", "xenium", "outlier_ids.csv"), row.names=FALSE)