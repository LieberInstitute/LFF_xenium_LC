library(data.table)
library(SpatialExperiment)
library(scuttle)
library(Banksy)
library(BiocParallel)
## faster UMAP in uwot 0.2 (used by banksy UMAP)
library(RcppHNSW)
library(here)

setwd("/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/")

# MAY NEED TO BE RERUN DEPENDING ON https://github.com/prabhakarlab/Banksy/issues/64

lcx <- readRDS("processed-data/04_Banksy_clusterings/01-staggeredCoord_spe_genetargOnly_NONlogNorm.RDS")

# get SLURM job #, i.e. row from the parameter table to use
job <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# get the parameter table
parms <- fread("processed-data/04_Banksy_clusterings/01-BanksyParamSets.txt")

# get parameter set for this job #
gabor <- parms[job,gabor]
kgeom <- parms[job,kgeom]
lam <- parms[job,lam]

# specify the non-log normalized counts assay to use; set a consistent seed
useassay <- "normcounts"
monocot = 42 # seed

# aaand run Banksy. this would be the same as NOT rearranging coords, then running on an SPE for each sample individually and joining them all back together after this step into one  SPE again. I think.

lcx <- computeBanksy(lcx, 
                      assay_name = useassay,
                      compute_agf = gabor,
                      k_geom = kgeom,
                      seed=monocot,
                      parallel = TRUE,
                      num_cores = 8)

## run Banksy PCA for subseequent input to clustering
## CRITICAL HERE: specify "group" as the sample_id -- i.e., each
## sample's Banksy PCA will be computed separately, making this
## equivalent to running Banksy on a list of SPEs with each entry as
## one sample, the other approach their vignettes uses
lcx <- runBanksyPCA(lcx,
                         assay_name=useassay,
                         use_agf=gabor,
                         lambda=lam,
                         group="sample_id",
                         seed=monocot)

## make sure to grab the BanksyPCA in case we want to subsequently Harmony it.
## we only need to save this once since we can run louvain at all desired
## resolutions on this. clusterBanksy does not parallelize so we run all resolutions
## in parallel instead.

pcaoutname <- paste0("PCA_M",ifelse(gabor,1,0),"_lam",lam)
pcaout <- reducedDims(lcx)[[pcaoutname]]
pcaoutpath <- paste0(
        "processed-data/04_Banksy_clusterings/01-banksy_clustruns_outs/banksyPCAs/",
        pcaoutname,"_kgeom",kgeom,".RDS")

if(!file.exists(pcaoutpath))
saveRDS(pcaout,paste0(
        "processed-data/04_Banksy_clusterings/01-banksy_clustruns_outs/banksyPCAs/",
        pcaoutname,"_kgeom",kgeom,".RDS"))

## run Banksy clustering -- louvain clustering (default method is leiden).
## several unwritten parameters will default to their defaults here. 
## and this will take a longgggg time, though louvain seems snappier.
## note that by default clusterBanksy would use k=50 for the 
## knn graph upstream of louv/leid which is a silly big number here.
## a two sample test run with k=50 at resolution 1 -> 19 clusters;
## k_neighb=25, res 1 -> 23 clusters, so similar # of clusters. we'll use
## k_neighbors=25, resolutions 0.5, 1, and 2.

### since louvain clustering itself cannot be parallelized and this job is sitting
### on 8 cpus, we can run all three resolutions at once. to avoid stupid collisions
### in memory, make copies of the SPE (its not THAT big) to ensure each process
### has its own copy. do this with base parallel:: , i don't trust biocp with
### SPE objects
lcx3res <- list(lcx,lcx,lcx)
library(parallel)
bclusts <- mcmapply(s=lcx3res,res=c(0.5,1,2),mc.cores=3,SIMPLIFY=FALSE,function(s,res){
        o <- Banksy::clusterBanksy(s,
              assay_name=useassay,
              use_agf = gabor,
              lambda = lam,
              k_neighbors = 25,
              resolution = res,
              group="sample_id",
              seed = monocot,algo = "louvain",)

        ## make the clustering result non-numeric for ease of plotting later;
        ## extract the clustering result, rather than saving a whole SPE for this one clustering result of manyy being tested.
        # the clustering result is a new column tacked onto colData.

        keepname <- c("rn",colnames(colData(o))[ncol(colData(o))])
        outtab <- as.data.table(colData(o),keep.rownames=T)[,..keepname]
        setnames(outtab,2,"clusts")
        outtab[,clusts:=paste0("X",as.character(clusts))]
        # append res (its already there, we fix it in the final table)

        setnames(outtab,2,paste0(keepname[2],"_kgeom",kgeom,"_res",res))
        
        # this column name will tell us each of the clustering parameters used for this run
        # e.g., clust_M0_lam0_k25_res2 means
        # 1. M=0 (no gabor), 
        # 2. lambda = 0,
        # 3. k25 means k neighbors for leiden clustering = 25, which we didn't vary here
        # 4. res2 means clustering resolution 2, 
        # 5. and then we add kgeom since we varied that parameter too.
        
return(outtab)})

# concat and save the cluster assignments
stopifnot(nrow(bclusts[[1]])==nrow(bclusts[[2]])&nrow(bclusts[[1]])==nrow(bclusts[[3]]))
outtab <- merge.data.table(bclusts[[1]],bclusts[[2]],by="rn")
outtab <- merge.data.table(outtab,bclusts[[3]],by="rn")

## get rid of the redundant, less visible resolution info
setnames(outtab,gsub(names(outtab),pattern="k25_res._",replacement="k25_"))

outtabname <- gsub(pcaoutname,pattern="PCA_",replacement="")
outtabname <- paste0(outtabname,"_kgeom",kgeom,"_louv05_louv1_louv2.txt")

fwrite(outtab,
       paste0("processed-data/04_Banksy_clusterings/01-banksy_clustruns_outs/",
              names(outtab)[2],
              ".txt"),
       sep='\t',
       quote=F,
       row.names=F,
       col.names=T)

## print reproducibility info
sessionInfo()
