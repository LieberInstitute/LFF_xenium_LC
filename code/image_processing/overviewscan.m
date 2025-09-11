Md= '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
Br = '/raw-data/xenium/output-XETG00558__0068654__Br6538__20250501__172909/aux_outputs/';
pix_json = jsondecode(fileread(fullfile(Md, Br, 'overview_scan_fov_locations.json')));
micron_json = jsondecode(fileread(fullfile(Md, Br, 'morphology_fov_locations.json')));

% Extract FOV locations
fovs = pix_json.fov_locations;
fovNames = fieldnames(fovs);


% Preallocate
xVals = zeros(length(fovNames), 1);
yVals = zeros(length(fovNames), 1);
widths = zeros(length(fovNames), 1);
heights = zeros(length(fovNames), 1);

% Loop through FOVs
for i = 1:length(fovNames)
    fov = fovs.(fovNames{i});
    xVals(i) = fov.x;
    yVals(i) = fov.y;
    widths(i) = fov.width;
    heights(i) = fov.height;
end

% Compute bounding box from x, y and their widths/heights
x_min = min(xVals);
x_max = max(xVals + widths);
y_min = min(yVals);
y_max = max(yVals + heights);

Oscan = imread(fullfile(Md, Br, 'overview_scan.png'));
cellmaskL = Oscan(y_min:y_max,x_min:x_max);
imwrite(cellmaskL,fullfile(Md,'processed-data/xenium_imageProcessing/overviewscan_Br6538.png'))

%%%% mpp extraction %%%

pix_dims = pix_json.fov_locations.AB10;      % in pixels
micron_dims = micron_json.fov_locations.AB10;  % in microns

% Compute microns per pixel
mpp_x = micron_dims.width / pix_dims.width;
mpp_y = micron_dims.height / pix_dims.height;

fprintf('Microns per pixel:\n  X: %.4f µm/pixel\n  Y: %.4f µm/pixel\n', ...
        mpp_x, mpp_y);

%Microns per pixel:
%  X: 7.1477 m/pixel
%  Y: 7.0566 m/pixel

[H_px, W_px, ~] = size(cellmaskL);

% Compute physical size in microns
width_um  = W_px * mpp_x;
height_um = H_px * mpp_y;

fprintf('Size of cellmaskL:\n  Width: %.2f µm (%d px)\n  Height: %.2f µm (%d px)\n', ...
        width_um, W_px, height_um, H_px);

%Size of cellmaskL:
%  Width: 6297.15 m (881 px)
%  Height: 2999.06 m (425 px)