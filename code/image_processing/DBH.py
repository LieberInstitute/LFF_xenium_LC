import scanpy as sc
import matplotlib.pyplot as plt

# Path to the H5 file
h5_file_new = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/code/xenium_imageProcessing/xenium_ranger/xeniumranger_NM_DAPI/outs/cell_feature_matrix.h5'
h5_file_old = '/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/_xenium_LIBD/250505_run-2/output-XETG00558__0068654__Br6297__20250501__172909/cell_feature_matrix.h5'

# Load the 10x H5 file into an AnnData object
adata_orig = sc.read_10x_h5(h5_file_old)
adata_new = sc.read_10x_h5(h5_file_new)

import pandas as pd

cellmeta_orig = pd.read_parquet('/dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/_xenium_LIBD/250505_run-2/output-XETG00558__0068654__Br6297__20250501__172909/cells.parquet')
cellmeta_new  = pd.read_csv('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/code/xenium_imageProcessing/xenium_ranger/xeniumranger_NM_DAPI/outs/cells.csv')

dbh_orig = adata_orig[:, 'DBH'].X.toarray().flatten()
dbh_new  = adata_new[:, 'DBH'].X.toarray().flatten()

total_dbh_orig = dbh_orig.sum()
total_dbh_new  = dbh_new.sum()

print("Original DBH+ cells:", (dbh_orig > 0).sum())
print("New DBH+ cells:", (dbh_new > 0).sum())
print("Original max area:", cellmeta_orig['nucleus_area'].max())
print("New max area:", cellmeta_new['nucleus_area'].max())
print(f"Original cell count: {adata_orig.n_obs}")
print(f"New cell count: {adata_new.n_obs}")
print("Total DBH counts assigned to cells (original):", total_dbh_orig)
print("Total DBH counts assigned to cells (new):", total_dbh_new)

if 'cell_id' in cellmeta_new.columns:
    cellmeta_new = cellmeta_new.set_index('cell_id')
    
if 'cell_id' in cellmeta_orig.columns:
    cellmeta_orig = cellmeta_orig.set_index('cell_id')

adata_new.obs = adata_new.obs.join(cellmeta_new, how='left')
adata_orig.obs = adata_orig.obs.join(cellmeta_orig, how='left')
  
plt.figure(figsize=(6,6))
plt.scatter(adata_orig.obs['x_centroid'], adata_orig.obs['y_centroid'], c=(dbh_orig > 0), cmap='Greys', s=2, label='Original')
plt.scatter(adata_new.obs['x_centroid'], adata_new.obs['y_centroid'], c=(dbh_new > 0), cmap='Oranges', s=2, label='New')
plt.legend()
plt.title('DBH+ Cell Locations')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()