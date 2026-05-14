baseHE = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing';
baseMask = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing_new';

brains = dir(fullfile(baseMask, 'Br*'));
brains = brains([brains.isdir]);

for i = 2:length(brains)

    brain = brains(i).name;
    fprintf('Processing %s\n', brain)

    he_file = fullfile(baseHE, brain, 'HE_registered_to_DAPI.tif');
    mask_file = fullfile(baseMask, brain, 'DAPImasks.mat');

    if ~isfile(he_file)
        fprintf('  Missing HE image: %s\n', he_file)
        continue
    end

    if ~isfile(mask_file)
        fprintf('  Missing DAPImasks.mat: %s\n', mask_file)
        continue
    end

    % Load HE
    HE = imread(he_file);

    % Load masks
    S = load(mask_file);
    cellmask = S.cellmask;
    nucmask  = S.nucmask;

    % Make binary masks from labeled masks
    cell_bin = cellmask > 0;
    nuc_bin  = nucmask > 0;

    % Boundaries
    cell_perim = bwperim(cell_bin);
    nuc_perim  = bwperim(nuc_bin);

    % Overlay
   	out = uint8(imoverlay(HE, cell_perim, [0 0 0]));
	out = uint8(imoverlay(out, nuc_perim, [0 1 0]));
   
    % Save
    out_file = fullfile(baseMask, brain, 'HE_DAPI_cell_nuc_overlay.png');
    imwrite(imresize(out, 0.5), out_file);

    fprintf('  Saved: %s\n', out_file)
end