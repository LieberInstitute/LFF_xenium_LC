D = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/raw-data/Images/'; 
myfiles = dir(fullfile(D, '*.tif'));
matchingFiles = myfiles(contains({myfiles.name}, {'_A1.tif', '_B1.tif', '_C1.tif', '_D1.tif'}));
O = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/Images/NMseg/';
for i = 1:numel(matchingFiles)

he = imread(fullfile(matchingFiles(i).folder, matchingFiles(i).name));
tic
disp('Performing rgb to Lab color space conversion')
lab_he = rgb2lab(he); % convert from rgb color space to Lab color space
toc

ab = lab_he(:,:,2:3); % extract a*b color space from Lab
ab = im2single(ab);
tic
disp('Applying Kmeans')
N=2;
pixel_labels = imsegkmeans(ab,N,'NumAttempts',3); % apply Kmeans
toc

parfor ic = 1:N
mask{ic} = pixel_labels==ic;
cluster{ic} = he .* uint8(mask{ic});
imwrite(imresize(cluster{ic}, 0.5), fullfile(O, [matchingFiles(i).name(1:end-4), num2str(ic),'.png']))
end

NM = cluster{2};
save(fullfile(O, [matchingFiles(i).name(1:end-4), '.mat']),'NM','-v7.3')

disp(['done ',matchingFiles(i).name])
end