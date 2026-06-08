library(data.table)
library(SpatialExperiment)
library(imager)

## >= and <= comparators that account for floating point error
greaterOrEqual <- function(x,y) {
    precision <- sqrt(.Machine$double.eps)
    (x >= y) | (abs(x-y) <= precision)
}

lessOrEqual <- function(x,y){
    precision <- sqrt(.Machine$double.eps)
    (x <= y) | (abs(x-y) <= precision)
}

library(parallelly)
library(parallel)

setDTthreads(1)

setwd("/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC")


## load filtered xenium data
spe <- readRDS("processed-data_newseg/02_QC/02b-spe_filtered_genetargonly_log-and-nonlog-norms.RDS")

## get image filepaths
brainz <- list.files("processed-data/xenium_imageProcessing_new/registrations/",pattern="_registered.png",recursive = FALSE)
brainz.dt <- as.data.table(cbind(brainz,gsub(brainz,pattern="-L",replacement="L")))
setnames(brainz.dt,c("dir","brnum"))
brainz.dt[,brnum:=gsub(brnum,pattern="^(Br.*)_registered\\.png",replacement="\\1")]
brainz.dt <- brainz.dt[brnum %in% spe$brnum]

stopifnot(all(unique(spe$brnum) %in% brainz.dt$brnum))
fwrite(brainz.dt,"processed-data_newseg/02_QC/05-tmpbrainzdt.txt",sep='\t',quote=F,row.names=F,col.names=T)

# we can't parallelize something png-y here. so iterate...
if(!dir.exists("processed-data/xenium_imageProcessing_new/HE_DAPIregist_lowres")){
  dir.create("processed-data/xenium_imageProcessing_new/HE_DAPIregist_lowres")  
}

b<-1

for (b in c(1:nrow(brainz.dt))){
    samp <- brainz.dt[b,brnum]
    fullres.path <- paste0("processed-data/xenium_imageProcessing_new/registrations/",brainz.dt[b,dir])
    
    ## make an output directory for our files
    if (!dir.exists(paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_lowres/",samp))){
      dir.create(paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_lowres/",samp))
    }
    
    outdir <- paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_lowres/",samp,"/")
    
    # load the full res HE-DAPI registered image
    fullres <- load.image(file.path(fullres.path))
    
    ## calculate scale factors for a 800 px max-dimension image
    lowres_max_size=800
    hires_scalef=0.2125# um/px

    low_over_hi <- lowres_max_size / max(dim(fullres)[seq(2)])
    lowres_scalefact <- low_over_hi/hires_scalef
    sr_json <- list(
        tissue_hires_scalef = hires_scalef,
        tissue_lowres_scalef = lowres_scalefact,
        spot_diameter_fullres = (5/hires_scalef)) # 5um nuclei / (.2125 um/px) = px 
    
    outimg <- imager::resize(
    fullres,
    as.integer(low_over_hi * dim(fullres)[1]),
    as.integer(low_over_hi * dim(fullres)[2]))
    
    imager::save.image(outimg,paste0(outdir,"tissue_lowres_image.png"))
    
    write(
        rjson::toJSON(sr_json),
        file.path(paste0(outdir,"scalefactors_json.json")))
    
    rm(fullres,outimg,sr_json,lowres_scalefact,low_over_hi,outdir,samp,dir)
    gc(full=T)
}

spes <- lapply(unique(spe$brnum),function(b){spe[,spe$brnum==b]})
names(spes) <- unique(spe$brnum)

newspes <- mcmapply(d=spes,n=names(spes),mc.cores=7,SIMPLIFY=FALSE,function(d,n){
    img <- readImgData(path=paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_lowres/",n,"/"),sample_id = n,scaleFactors = paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_lowres/",n,"/scalefactors_json.json"))
    
    newspe <- SpatialExperiment(assays=list(counts=counts(d)),sample_id = n,spatialCoords = spatialCoords(d),imgData = img,colData=colData(d),rowData=rowData(d))

    return(newspe)    
})
names(newspes) <- names(spes)

i<-1
for (i in c(1:length(newspes))){
    
    if(i==1){
        outspe <- newspes[[i]]
    }
    else {
        outspe <- cbind(outspe,newspes[[i]])
    }
}
stopifnot(dim(outspe)==dim(spe))
dim(outspe)
length(unique(outspe$brnum))

#save
saveRDS(outspe,"processed-data_newseg/02_QC/05a-spe_filtered_genetargonly_log-and-nonlog-norms_withLowresHEs.RDS")

# seshinf
sessionInfo()
sessioninfo::session_info()
