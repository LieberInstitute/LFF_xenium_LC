import pandas as pd

df = pd.read_parquet("/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/raw-data/xenium/output-XETG00558__0068654__Br6538__20250501__172909/transcripts.parquet")
df.to_csv("transcripts.csv", index=False)  # or df.to_hdf or df.to_pickle

genes_of_interest = ['TH', 'DBH', 'PHOX2A', 'PHOX2B', 'SLC6A2']
df_subset = df[df['feature_name'].isin(genes_of_interest)]
df_subset.to_csv("/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/Br6538_transcripts_subset.csv", index=False)