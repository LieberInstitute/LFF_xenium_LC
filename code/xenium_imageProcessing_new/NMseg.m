function NMseg(brnum)
% K-means (Lab a*b) clustering of a registered HE image.
% Writes: Kmeans_labels.tif, Kmeans_masks.mat, Kmeans_cluster.mat, Kmeans_cluster.png

    fprintf('Brnum: %s\n', brnum);

    baseDir = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing_new/registrations';
    outDir  = fullfile(baseDir, brnum);
    inPath  = fullfile(baseDir, [brnum, '_HE.png']);

    if ~exist(inPath, 'file')
        error('Input not found: %s', inPath);
    end
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Load image
    he = imread(inPath);           % likely uint16 RGB
    if size(he,3) ~= 3
        error('Expected RGB image (MxNx3), got %s', mat2str(size(he)));
    end

    % Convert to Lab (ab) for color clustering
    tic; disp('Converting RGB -> Lab...'); 
    lab_he = rgb2lab(he);          % uint16 OK; rgb2lab handles it
    ab = lab_he(:,:,2:3);
    ab = im2single(ab);            % imsegkmeans expects floating
    toc

    % If very large, consider downsampling before k-means (optional)
    % scale = 1; if max(size(he,1), size(he,2)) > 12000, scale = 0.5; end
    % if scale ~= 1, ab_small = imresize(ab, scale, 'bilinear'); else, ab_small = ab; end

    % K-means clustering
    N = 2;
    tic; disp('Running k-means on a*b...'); 
    pixel_labels = imsegkmeans(ab, N, 'NumAttempts', 3);
    toc

    % Prepare outputs
    mask    = cell(1,N);
    cluster = cell(1,N);

    heClass = class(he);  % remember original class for clean masking
    for ic = 1:N
        mask{ic} = (pixel_labels == ic);
        % cast mask to the image class before multiplication
        switch heClass
            case 'uint16', mskTyped = uint16(mask{ic});
            case 'uint8',  mskTyped = uint8(mask{ic});
            otherwise,     mskTyped = cast(mask{ic}, heClass);
        end
        cluster{ic} = he .* mskTyped;
    end

    save(fullfile(outDir, 'Kmeans_cluster.mat'), 'cluster', '-v7.3');

    % Save a quicklook PNG (convert safely to 8-bit)
    % Use im2uint8 so the dynamic range is preserved sensibly
    C1 = im2uint8(rescale(cluster{1}));
    C2 = im2uint8(rescale(cluster{2}));
    imwrite(imresize([C1, C2], 0.2), fullfile(outDir, 'Kmeans_cluster.png'));

    disp('Done.');
end