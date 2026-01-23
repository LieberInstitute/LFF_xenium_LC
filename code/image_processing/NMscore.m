cd '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC'

	samples = readtable('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code/image_processing/inputs.txt')
%% samples
files = dir(fullfile(pwd, '/processed-data/xenium_imageProcessing/'));
myfiles = files(cellfun(@(x) length(x) == 17, {files.name}));

%% bernie annotation based BG 
for i= 1:numel(myfiles)
fname = myfiles(i).name(1:end-4);
disp(fname);
img = imread([pwd, '/raw-data/Images/',fname,'.tif']);

NMseg_dir = fullfile(pwd, '/processed-data/Images/NMseg/');
load([NMseg_dir, fname, 'NMseg_clean.mat'])

img1 = mat2gray(rgb2gray(img));
BW = img1<0.8&img1>0.1;     
img1(~BW)=0; 
	
BG_mask = img1; BG_mask(NM) = 0; 
HE = mean(BG_mask>0);

NM_mask = img1; NM_mask(~NM) = 0; 
nm = NM_mask(NM_mask>0);
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
SNM_score1 = Pi/Mt;


