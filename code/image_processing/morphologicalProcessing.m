Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
od = '/processed-data/xenium_imageProcessing/';
brain = 'Br6297';
disp(brain);

load(fullfile(Md,od,brain, 'NMseg.mat'));
BW1 = zeros(size(BW));
BW1(8000:21000,2000:12000) = BW(8000:21000,2000:12000);

cleaned = imopen(BW1, strel('disk', 3)); 
cleanedN = imclose(cleaned, strel('disk', 3)); 
stats = regionprops(logical(cleanedN), 'Area'); 
%keptSize = ([stats.Area] >= 80); % & ([stats.Area] <= 5000);
%kept = find(keptSize);   
%NM = ismember(labelmatrix(bwconncomp(logical(cleanedN))), kept);
NM = bwareaopen(cleanedN, 500);
NM(20000:21000,2000:4000) = 0;
%img = im2double(imread(fullfile(Md,od,brain,[brain,'_HE_aligned.png'])));
img_overlay = img;
nm_mask = NM > 0;

% Set yellow color: R and G to 1, B to 0 where nm_mask is true
img_overlay(:,:,1) = img(:,:,1) .* ~nm_mask + nm_mask * 1;  % Red channel
img_overlay(:,:,2) = img(:,:,2) .* ~nm_mask + nm_mask * 1;  % Green channel
img_overlay(:,:,3) = img(:,:,3) .* ~nm_mask + nm_mask * 0;  % Blue channel

% Display and save
figure, imshow(img_overlay(8000:21000,2000:12000,:))
save(fullfile(Md,od,brain, 'NMseg_clean.mat'), 'NM', '-v7.3');



