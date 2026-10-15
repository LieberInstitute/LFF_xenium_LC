library(data.table)
library(SpatialExperiment)
library(qs2)
library(parallel)
library(parallelly)

setDTthreads(1)

setwd("/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC")

spe <- readRDS("processed-data_newseg/02_QC/02b-spe_filtered_genetargonly_log-and-nonlog-norms.RDS")

spes <- lapply(unique(spe$brnum),function(b){spe[,spe$brnum==b]})
names(spes) <- unique(spe$brnum)

newspes <- mcmapply(d=spes,n=names(spes),mc.cores=4,SIMPLIFY=FALSE,function(d,n){
    img <- readImgData(path=paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_midres/",n,"/"),sample_id = n,scaleFactors = paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_midres/",n,"/scalefactors_json.json"))

    newspe <- SpatialExperiment(assays=list(counts=counts(d)),sample_id = n,spatialCoords = spatialCoords(d),imgData = img,colData=colData(d),rowData=rowData(d))
})

names(newspes) <- names(spes)

### save single-sample SPEs with larger H&E image for Xenium viewer
mapply(d=spes,n=names(spes),SIMPLIFY=FALSE,function(d,n){
    qs_save(d,paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_midres/",n,".qs2"))
    # return(newspe)
})

# seshinf
sessionInfo()
sessioninfo::session_info()
