import pandas as pd
import glob
import os

base_dir = "/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing"

# find all Br**** folders
br_folders = glob.glob(os.path.join(base_dir, "Br*"))

genes_of_interest = ['DBH', 'P2RY12', 'FGFR', 'MOG', 'FLT1', 'FCGR1A', 'P2RY13', 'MAG', 'APOLD1']

for br_path in br_folders:
    br_name = os.path.basename(br_path)  # e.g. Br1234
    parquet_path = os.path.join(
        br_path,
        f"xeniumranger_NM_DAPI_{br_name}",
        "outs",
        "transcripts.parquet"
    )
    if not os.path.exists(parquet_path):
        print(f"Skipping {br_name} (no parquet found)")
        continue
    print(f"Processing {br_name}")
    df = pd.read_parquet(parquet_path)
    # save full table (optional)
    #df.to_csv(os.path.join(br_path, "transcripts.csv"), index=False)
    # subset genes
    df_subset = df[df['feature_name'].isin(genes_of_interest)]
    
    df_subset.to_csv(
        os.path.join(br_path, "transcripts_subset.csv"),
        index=False
)

#### 