cd '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC'

	fname = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code/image_processing/inputs.txt';

	opts = detectImportOptions(fname, 'Delimiter', '\t');
	opts.DataLines = [1 Inf];     % <- don’t skip first 4 lines
	opts.CommentStyle = {};       % <- don’t drop lines starting with #, etc.
	samples = readtable(fname, opts);

for i= 1:size(samples,1)	
	disp(samples.Var1{i})
		
img = imread(fullfile(pwd, '/processed-data/xenium_imageProcessing/', samples.Var1{i}, 'HE_registered_to_DAPI.tif'));
load(fullfile(pwd, '/processed-data/xenium_imageProcessing/', samples.Var1{i}, 'NMnofolds.mat'))

img1 = mat2gray(rgb2gray(img));
BW = img1<0.8&img1>0.1;     
img1(~BW)=0; 
	
BG_mask = img1; BG_mask(BW_nofolds) = 0; temp=BG_mask(BG_mask>0);
HE = mean(BG_mask(:));

NM_mask = img1; NM_mask(~BW_nofolds) = 0; 
nm = NM_mask(NM_mask>0);
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
SNM_score1 = Pi/Mt;


