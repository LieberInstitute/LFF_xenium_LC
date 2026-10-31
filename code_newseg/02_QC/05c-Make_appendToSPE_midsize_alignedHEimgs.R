library(data.table)
library(SpatialExperiment)
library(imager)
library(qs2)

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
if(!dir.exists("processed-data/xenium_imageProcessing_new/HE_DAPIregist_midres")){
  dir.create("processed-data/xenium_imageProcessing_new/HE_DAPIregist_midres")  
}

b<-1

for (b in c(1:nrow(brainz.dt))){
    samp <- brainz.dt[b,brnum]
    fullres.path <- paste0("processed-data/xenium_imageProcessing_new/registrations/",brainz.dt[b,dir])
    
    ## make an output directory for our files
    if (!dir.exists(paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_midres/",samp))){
      dir.create(paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_midres/",samp))
    }
    
    outdir <- paste0("processed-data/xenium_imageProcessing_new/HE_DAPIregist_midres/",samp,"/")
    
    # load the full res HE-DAPI registered image
    fullres <- load.image(file.path(fullres.path))
    
    ## calculate scale factors for a 3500 px max-dimension image
    lowres_max_size=3500
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
    # we have to save this as the "lowres" image for consistency with what
    # spatialexperiment looks for when loading image data.
    imager::save.image(outimg,paste0(outdir,"tissue_lowres_image.png"))
    
    write(
        rjson::toJSON(sr_json),
        file.path(paste0(outdir,"scalefactors_json.json")))
    
    rm(fullres,outimg,sr_json,lowres_scalefact,low_over_hi,outdir,samp,dir)
    gc(full=T)
}

# seshinf
sessionInfo()
sessioninfo::session_info()
