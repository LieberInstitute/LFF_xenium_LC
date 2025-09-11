Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
od = '/processed-data/xenium_imageProcessing/';
brain = 'Br6297';
load(fullfile(Md, od, brain, 'nuc.mat'))
load(fullfile(Md, od, brain, 'NMseg_clean.mat'))
NM_labeled = bwlabel(NM);
%% seperate NM labels from nucmask labels
%Offset NM labels to avoid overlap with nucmask labels
max_nuc_label = max(nucmask(:));
NM_labeled(NM_labeled > 0) = NM_labeled(NM_labeled > 0) + max_nuc_label;

combined_mask = nucmask;
combined_mask(NM_labeled > 0) = NM_labeled(NM_labeled > 0);

tiff_filename = fullfile(Md, od, brain, 'combined_nucmask.tif');

% Create a Tiff object
t = Tiff(tiff_filename, 'w');

% Set TIFF tags for 32-bit image
tagstruct.ImageLength = size(combined_mask, 1);
tagstruct.ImageWidth = size(combined_mask, 2);
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Compression = Tiff.Compression.None;
t.setTag(tagstruct);
t.write(uint32(combined_mask));
t.close();


imwrite(uint32(combined_mask), tiff_filename, 'Compression', 'none');

%% save as Numpy array
combined_mask_uint32 = uint32(combined_mask);
np_array = py.numpy.array(combined_mask_uint32);

output_file = fullfile(Md, od, brain, 'combined_nucmask.npy');
py.numpy.save(output_file, np_array);


%% combine NM regions with nucmask only if there is 30% overlap
NM_labeled = bwlabel(NM);
combined_mask = nucmask;
max_label = max(nucmask(:));

NM_stats = regionprops(NM_labeled, 'PixelIdxList');
nucmask_stats = regionprops(nucmask, 'PixelIdxList');

for i = 1:numel(NM_stats)
    nm_pixels = NM_stats(i).PixelIdxList;
    overlapping_labels = nucmask(nm_pixels);
    overlapping_labels(overlapping_labels == 0) = [];
    
    if isempty(overlapping_labels)
        % No overlap, assign new label
        max_label = max_label + 1;
        combined_mask(nm_pixels) = max_label;
        continue;
    end

    % Count how many pixels overlap with each nucmask label
    [ulabels, ~, idx] = unique(overlapping_labels);
    counts = accumarray(idx, 1);
    
    assigned = false;
    for j = 1:numel(ulabels)
        nuc_label = ulabels(j);
        overlap_count = counts(j);
        nuc_size = numel(nucmask_stats(nuc_label).PixelIdxList);
        
        overlap_ratio = overlap_count / nuc_size;
        if overlap_ratio >= 0.3
            combined_mask(nm_pixels) = nuc_label;
            assigned = true;
            break;
        end
    end

    if ~assigned
        max_label = max_label + 1;
        combined_mask(nm_pixels) = max_label;
    end
end

rgb_nucmask = label2rgb(uint16(nucmask), 'jet', 'k', 'shuffle');
rgb_NMnucmask = label2rgb(uint16(combined_mask), 'jet', 'k', 'shuffle');

