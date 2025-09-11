Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
od = '/processed-data/xenium_imageProcessing/';
brain = 'Br6297';
filename = fullfile(Md, od, brain, 'xeniumranger_NM_DAPI', 'outs', 'transcripts.parquet');
T = parquetread(filename);