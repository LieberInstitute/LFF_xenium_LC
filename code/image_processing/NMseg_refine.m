Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
od = '/processed-data/xenium_imageProcessing/';
brain = 'Br6297';
fname = 'Br6297_Kmeans.mat';
disp(fname);

   load(fullfile(Md,od,brain,fname));
   he = imread(fullfile(Md,od,brain,[brain,'_HE_aligned.png']));
   NM = cluster{1};
   Ie = rgb2gray(NM);
   Ie = im2double(Ie);
   BW = zeros(size(Ie));
   BW(Ie>0 & Ie<0.25) = 1;
   
   save(fullfile(Md,od,brain, 'NMseg.mat'), 'BW', '-v7.3');

   BW = 255*cat(3, BW, BW, BW);
   bw = BW(9000:11000, 8000:10000,:);
   nm = NM(9000:11000, 8000:10000,:);
   img = he(9000:11000, 8000:10000,:);

   imwrite([img,nm,bw], fullfile(Md,od,brain, 'NMseg.png'))
