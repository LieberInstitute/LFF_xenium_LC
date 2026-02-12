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
res <- parms[job,res]

# specify the non-log normalized counts assay to use; set a consistent seed
useassay <- "normcounts"
monocot = 42 # seed

# aaand run Banksy. this would be the same as NOT rearranging coords, then running on an SPE for each sample individually and joining them all back together after this step into one  SPE again. I think.

lcx <- computeBanksy(lcx, 
                      assay_name = useassay,
                      compute_agf = gabor,
                      k_geom = kgeom,
                      seed=monocot)

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

## run Banksy clustering -- Leiden clustering (default method).
## several unwritten parameters will default to their defaults here. 
## and this will take a longgggg time.
lcx <- clusterBanksy(lcx,
                      assay_name=useassay,
                      use_agf = gabor,
                      lambda = lam,
                      resolution = res,
                      group="sample_id",
                      seed = monocot)

## make the clustering result non-numeric for ease of plotting later;
## extract the clustering result, rather than saving a whole SPE for this one clustering result of manyy being tested.
# the clustering result is a new column tacked onto colData.

keepname <- c("rn",colnames(colData(lcx))[ncol(colData(lcx))])
outtab <- as.data.table(colData(lcx),keep.rownames=T)[,..keepname]
setnames(outtab,2,"clusts")
outtab[,clusts:=paste0("X",as.character(clusts))]
setnames(outtab,2,paste0(keepname[2],"_kgeom",kgeom)) # restore the previous name, with the kgeom parameter added.

# this column name will tell us each of the clustering parameters used for this run
# e.g., clust_M0_lam0_k50_res2 means
# 1. M=0 (no gabor), 
# 2. lambda = 0,
# 3. k50 means k neighbors for leiden clustering = 50, which we didn't vary here
# 4. res2 means clustering resolution 2, 
# 5. and then we add kgeom since we varied that parameter too.

# save the cluster assignments
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
