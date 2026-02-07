setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC')

# Load libraries
library(data.table)
library(RcppML)
library(SpatialExperiment)
library(scuttle)
library(Matrix)
library(singlet)

# Load Visium data and annots
lcv <- readRDS("../LFF_spatial_LC/processed-data/06-countsOnly_noImgData_QCed_SPE_split_to_tissSections.RDS")

vdom <- readRDS("../LFF_spatial_LC/processed-data/07_featureSelection_dimred_harmony_clustering/04-Initial_louvLeid_res0.5-1-2_clusterings.RDS")
vdom <- vdom$HARMONYlmbna_HDG_SVG_2575
vdom <- vdom[,.(rn,snnHARMONYlmbna_HDG_SVG_2575_louv_res1)]
setnames(vdom,2,"clusid")

vanno <- fread("../LFF_spatial_LC/processed-data/08_validitycheck_25hdg75svg_louv1/10-25hdg75svg_louv1_annots.txt")
vdom <- merge.data.table(vdom,vanno[,.(clusid,anno)],by="clusid")
vdom <- DataFrame(vdom,row.names=vdom$rn)[colnames(lcv),]

colLabels(lcv) <- vdom$anno

rm(vdom,vanno)

unique(lcv$brnum)
lcv$brnum[lcv$brnum=="Br6119(re-dis)"] <- "Br6119"

## drop donors with limited LC and br 1691 (tentatively dropping from Xenium)
lcv <- lcv[,!(lcv$brnum %in% paste0("Br",c(1691,5517,5276,5712)))]


### we need logcounts on here to perform NMF
lcv <- scuttle::logNormCounts(lcv)

## run it with a lower tol 
options(RcppML.threads=15)
options(RcppML.verbose=TRUE)
setDTthreads(1,restore_after_fork=FALSE) # prevent resource competition
crossval <- RcppML::crossValidate(logcounts(lcv),k=seq(50,150,10),tol=1e-5,threads=7,alpha=0,verbose=T)


saveRDS(crossval, "processed-data/00_prelim_procdat_BJM/05_LabelTransfer/01a-Crossval_VisNMF_50to150.RDS")

library(ggplot2)
png("plots/00_prelimplotsBJM/05_Visdomaintransfer/01a-VisNMF_crossval_k50to150.png", height=4,width=8, unit="in",res=300)
ggplot(crossval,aes(x=k,y=value))+geom_point()+scale_x_continuous(breaks=seq(50,150,20))
dev.off()

sessionInfo()
