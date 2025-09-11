Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
od = '/processed-data/xenium_imageProcessing/';
brain = 'Br6297';

he = imread(fullfile(Md,od,brain,[brain,'_HE_aligned.png']));
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
imwrite(imresize(cluster{ic}, 0.5), fullfile(Md,od, [brain, '_Kmeans', num2str(ic),'.png']))
end

save(fullfile(Md,od, [brain, '_Kmeans', '.mat']),'cluster','-v7.3')

disp(['done ',brain])
