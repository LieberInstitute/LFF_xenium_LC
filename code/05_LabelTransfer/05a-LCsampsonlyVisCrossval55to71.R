setwd('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC')

# Load libraries
library(data.table)
library(RcppML)
library(SpatialExperiment)
library(scuttle)
library(Matrix)
library(singlet)

# Load Visium data and annots
lcv <- readRDS("../LFF_spatial_LC/processed-data/06-LCposDonorsOnly_countsOnly_noImgData_QCed_SPE_split_to_tissSections.RDS")

vdom <- fread("../LFF_spatial_LC/processed-data/07_2_LCsampsonly_featureSelection_dimred_harmony_clustering/04_2-LConly_HDG_SVG_2575_Louv1_clustering.txt")
setnames(vdom,2,"clusid")

vanno <- fread("../LFF_spatial_LC/processed-data/07_2_LCsampsonly_featureSelection_dimred_harmony_clustering/08_2-25hdg75svg_louv1_annots.txt")
vdom <- merge.data.table(vdom,vanno[,.(clusid,anno)],by="clusid")
vdom <- DataFrame(vdom,row.names=vdom$rn)[colnames(lcv),]

colLabels(lcv) <- vdom$anno

rm(vdom,vanno)

unique(lcv$brnum)
lcv$brnum[lcv$brnum=="Br6119(re-dis)"] <- "Br6119"

## drop donors with limited LC and br 1691 (tentatively dropping from Xenium)
lcv <- lcv[,!(lcv$brnum %in% paste0("Br",c(1691,5517,5276,5712)))]


### we need logcounts on here to perform NMF
lcv <- scater::computeLibraryFactors(lcv)
lcv <- scater::logNormCounts(lcv)

## run it with a lower tol 
options(RcppML.threads=15)
options(RcppML.verbose=TRUE)
setDTthreads(1,restore_after_fork=FALSE) # prevent resource competition
crossval <- RcppML::crossValidate(logcounts(lcv),k=seq(55,71,1),tol=1e-5,threads=15,alpha=0,verbose=T)


saveRDS(crossval, "processed-data/05_LabelTransfer/05a-Crossval_LCposOnly_Vis_NMF_55to71.RDS")

library(ggplot2)
png("plots/00_prelimplotsBJM/05_Visdomaintransfer/05a-LCposOnly_Vis_NMF_crossval_k55to71.png", height=4,width=8, unit="in",res=300)
ggplot(crossval,aes(x=k,y=value))+geom_point()+scale_x_continuous(breaks=seq(54,72,3))
dev.off()

sessionInfo()
