addpath('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/code/xenium_imageProcessing/bfmatlab')    

Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
img = imread(fullfile(Md,'/processed-data/xenium_imageProcessing/Br6538_HE.png'));
out_file = fullfile(Md,'processed-data/xenium_imageProcessing/Br6538_HE.ome.tif');

% Convert image to uint8 if needed
if ~isa(img, 'uint8')
    img = im2uint8(img);
end

% Write OME-TIFF
bfWriteImage(out_file, img, 'BigTiff', true);