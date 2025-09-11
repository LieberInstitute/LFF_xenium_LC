Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
od = '/processed-data/xenium_imageProcessing/';
brain = 'Br6297';

img_registered = imread(fullfile(Md,od,brain,[brain,'_HE_aligned.png']));
img = im2double(img_registered);

load(fullfile(Md, od,brain,'cell.mat'))
%scale = 0.84;
%cellmaskL = imresize(cellmask, scale, 'nearest');
cellmaskL=cellmask;

% Generate boundary mask
boundary_mask = bwperim(cellmaskL > 0);
thick_boundary = imdilate(boundary_mask, strel('disk', 2));

img(repmat(thick_boundary, [1 1 3])) = 0;

%% to plot boundaries of attached cells
% Round label image to ensure integers
cellmaskR = round(cellmaskL);

% Use dilation to find borders between different labeled regions
se = strel('disk', 1);
dilated = imdilate(cellmaskR, se);
boundary_mask = (dilated ~= cellmaskR) & (cellmaskR > 0);
thick_boundary = imdilate(boundary_mask, strel('disk', 2));

img(repmat(thick_boundary, [1 1 3])) = 0;
%imshow(img)
imwrite(im2uint8(img), fullfile(Md, od, brain, [brain, '_HE_DAPI_overlay.png']));

%% overlay NMseg
load(fullfile(Md,od,brain,'NMseg.mat'))
nm_mask = BW > 0;
img_overlay = img;

% Set yellow color: R and G to 1, B to 0 where nm_mask is true
img_overlay(:,:,1) = img(:,:,1) .* ~nm_mask + nm_mask * 1;  % Red channel
img_overlay(:,:,2) = img(:,:,2) .* ~nm_mask + nm_mask * 1;  % Green channel
img_overlay(:,:,3) = img(:,:,3) .* ~nm_mask + nm_mask * 0;  % Blue channel

% Optional: blend instead of overwrite
% alpha = 0.5;
% img_overlay(:,:,1) = img(:,:,1) .* (1 - alpha * nm_mask) + alpha * nm_mask;
% img_overlay(:,:,2) = img(:,:,2) .* (1 - alpha * nm_mask) + alpha * nm_mask;
% img_overlay(:,:,3) = img(:,:,3) .* (1 - alpha * nm_mask);  % stays dark in blue

% Display and save
imshow(img_overlay(9000:10000,8000:10000,:))
imwrite(im2uint8(img_overlay), fullfile(Md, od, brain, [brain, '_HE_DAPI_NM_overlay.png']));

%% iPSC grant overlay transcripts

T = readtable(fullfile(Md,od,'Br6538_transcripts_subset.csv'));
microns_per_pixel = 0.253;
T.x_pixel = round(T.x_location / microns_per_pixel);
T.y_pixel = round(T.y_location / microns_per_pixel);

% Subset to TH only
idx_TH = strcmp(T.feature_name, 'TH');

% Display and overlay
figure;
imshow(img);
hold on;
scatter(T.x_pixel(idx_TH), T.y_pixel(idx_TH), 0.03, 'r', 'filled');  % Red dots
title('TH Transcripts Overlay');
hold off;

%%%%% all transcripts
% Define target genes and distinct colors
target_genes = {'TH', 'DBH', 'SLC6A2'};
colors = [
    1 0 0;    % Red
    0 1 0;    % Green
    1 1 0     % Yellow
];

% Plot image
figure;
imshow(img);
hold on;

% Plot each gene
for i = 1:length(target_genes)
    gene = target_genes{i};
    gene_idx = strcmp(T.feature_name, gene);
    scatter(T.x_pixel(gene_idx), T.y_pixel(gene_idx), 0.01, colors(i,:), 'filled', ...
        'DisplayName', gene);
end

legend('Location', 'bestoutside');
title('Overlay of TH, DBH, SLC6A2 transcripts');
hold off;
