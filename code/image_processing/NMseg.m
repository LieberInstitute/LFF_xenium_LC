function NMseg(brnum)

disp(brnum)
Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/'
fname = fullfile(Md,brnum,'HE_registered_to_DAPI.tif')	
	
he = imread(fname);
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

disp('saving ...')
	
parfor ic = 1:N
mask{ic} = pixel_labels==ic;
cluster{ic} = he .* uint8(mask{ic});
imwrite(cluster{ic}, fullfile(Md, brnum, ['Kmeans_cluster', num2str(ic),'.png']))
end

end