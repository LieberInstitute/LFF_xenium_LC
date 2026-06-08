baseNew = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing_new';
baseSplit = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/split_samples';
outDir = fullfile(baseNew, 'registrations');

inputFile = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code/xenium_imageProcessing_new/inputs.txt';

C = readcell(inputFile, ...
    'FileType', 'text', ...
    'Delimiter', '\t');

T = cell2table(C);

for i = 2:height(T)

    brain = string(T.C1{i});
    slideDir = string(T.C2{i});
    heName = string(T.C3{i});

    fprintf('Processing %s\n', brain)

    maskFile = fullfile(baseNew, brain, 'DAPImasks.mat');
    heFile = fullfile(baseSplit, slideDir, heName);

    if ~isfile(maskFile)
        fprintf('  Missing mask: %s\n', maskFile)
        continue
    end

    if ~isfile(heFile)
        fprintf('  Missing HE: %s\n', heFile)
        continue
    end

    % Load nucmask
    load(maskFile, 'nucmask');

    % Make nucmask boundary
    nuc_bin = nucmask > 0;
    nuc_perim = bwperim(nuc_bin);

    % Save nucmask boundary PNG
    nucOut = fullfile(outDir, brain + "_nucmaskPerim.png");
    imwrite(uint8(nuc_perim) * 255, nucOut);

    % Load original HE
    HE = imread(heFile);

    % Save original HE as Br****_HE.png
    heOut = fullfile(outDir, brain + "_HE.png");

    if ~isa(HE, 'uint8')
        HE = im2uint8(HE);
    end

    imwrite(HE, heOut);

    fprintf('  Saved: %s\n', nucOut)
    fprintf('  Saved: %s\n', heOut)
end