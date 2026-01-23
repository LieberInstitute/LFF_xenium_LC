cd '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC'
%% load metadata and BB
load(fullfile(pwd, '/processed-data/Images/NMseg/Mdata.mat'))

%% load bernies annotations
annot = fullfile(pwd, '/processed-data/Images/06c-tissueOverlapping_spots_with_fullres-pixel-row-col-coords_and_section_annotations.txt');

tb = readtable(annot);
naIdx = strcmp(tb.section, 'NA');
naIndices = find(naIdx);
if size(naIndices, 1)> 0 , disp([naIndices, 'donot have annotations']); end

%% samples
files = dir(fullfile(pwd, '/raw-data/Images/*1.tif'));
myfiles = files(cellfun(@(x) length(x) == 17, {files.name}));

%% bernie annotation based BG 
for i= 1:numel(myfiles)
fname = myfiles(i).name(1:end-4);
disp(fname);
img = imread([pwd, '/raw-data/Images/',fname,'.tif']);

NMseg_dir = fullfile(pwd, '/processed-data/Images/NMseg/');
load([NMseg_dir, fname, 'NMseg_clean.mat'])

img1 = mat2gray(rgb2gray(img));
BW = img1 < 0.9;     

% bernie annotation based NMScore
df = tb(strcmp(tb.sample_id,fname) & strcmp(tb.section, 'section_1'), :);
x = df.pxl_col_in_fullres; y = df.pxl_row_in_fullres;
k = convhull(y,x);
Bmask = poly2mask(x(k), y(k), size(BW, 1), size(BW, 2)); % Black out everything outside the polygon
BmaskImg = bsxfun(@times, BW, cast(Bmask, class(BW)));
img1(~BmaskImg) = 0;

BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;

HE = mean(BG_mask(:));
nm = NM_mask(NM_mask>0);
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
SNM_score1 = Pi/Mt;

Mdata.NMscoreS1(Mdata.sample_id == fname) = SNM_score1;

img1 = mat2gray(rgb2gray(img));
df = tb(strcmp(tb.sample_id,fname) & strcmp(tb.section, 'section_2'), :);
x = df.pxl_col_in_fullres; y = df.pxl_row_in_fullres;
k = convhull(y,x);
Bmask = poly2mask(x(k), y(k), size(BW, 1), size(BW, 2)); % Black out everything outside the polygon
BmaskImg = bsxfun(@times, BW, cast(Bmask, class(BW)));
img1(~BmaskImg) = 0;

BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;

HE = mean(BG_mask(:));
nm = NM_mask(NM_mask>0);
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
SNM_score2 = Pi/Mt;
Mdata.NMscoreS2(Mdata.sample_id == fname) = SNM_score2;


% BB annotation based NMScore
img1 = mat2gray(rgb2gray(img));
roi = Mdata.BB1{Mdata.sample_id == fname};  
x = roi(1); y = roi(2); w = roi(3); h = roi(4);
mask = false(size(img,1), size(img,2));
mask(y:y+h-1, x:x+w-1) = true;  % Set ROI region to true
maskImg = bsxfun(@times, BW, cast(mask, class(BW)));
img1(~maskImg) = 0;

BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;

HE = mean(BG_mask(:));
nm = NM_mask(NM_mask>0);
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
NM_score = Pi/Mt;

Mdata.NMscoreB1(Mdata.sample_id == fname) = NM_score;
Mdata.NMscorebB1(Mdata.sample_id == fname) = mean(nm);

%bb2
img1 = mat2gray(rgb2gray(img));
roi = Mdata.BB2{Mdata.sample_id == fname};  
x = roi(1); y = roi(2); w = roi(3); h = roi(4);
mask = false(size(img,1), size(img,2));
mask(y:y+h-1, x:x+w-1) = true;  % Set ROI region to true
maskImg = bsxfun(@times, BW, cast(mask, class(BW)));
img1(~maskImg) = 0;

BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;

HE = mean(BG_mask(:));
nm = NM_mask(NM_mask>0);
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
NM_score = Pi/Mt;

Mdata.NMscoreB2(Mdata.sample_id == fname) = NM_score;
Mdata.NMscorebB2(Mdata.sample_id == fname) = mean(nm);

disp(i)
end

writetable(Mdata, fullfile(pwd,'processed-data/Images/NMscore.csv'));





%%%%%%%
% compare mean nm intensities before and after BG subtraction for bernie
% annotated regions

for i= 1:numel(myfiles)
fname = myfiles(i).name(1:end-4);
disp(fname);
img = imread([pwd, '/raw-data/Images/',fname,'.tif']);

NMseg_dir = fullfile(pwd, '/processed-data/Images/NMseg/');
load([NMseg_dir, fname, 'NMseg_clean.mat'])

img1 = mat2gray(rgb2gray(img));
BW = img1 < 0.9;     

% bernie annotation based NMScore
df = tb(strcmp(tb.sample_id,fname) & strcmp(tb.section, 'section_1'), :);
x = df.pxl_col_in_fullres; y = df.pxl_row_in_fullres;
k = convhull(y,x);
Bmask = poly2mask(x(k), y(k), size(BW, 1), size(BW, 2)); % Black out everything outside the polygon
BmaskImg = bsxfun(@times, BW, cast(Bmask, class(BW)));
img1(~BmaskImg) = 0;

BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;
nm = NM_mask(NM_mask>0);

% Pi = sum(nm);
% maxNM = max(nm);
% Mt = size(nm,1)*maxNM;
% NMscore1 = Pi/Mt;
% 
% Mdata.NMscore1(Mdata.sample_id == fname) = NMscore1;

HE = mean(BG_mask(:));
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
BGNMscore1 = Pi/Mt;
temp1 = sum(nm-HE < 0);

%Mdata.BGNMscore1(Mdata.sample_id == fname) = BGNMscore1;
Mdata.tempsize1(Mdata.sample_id == fname) = temp1;

img1 = mat2gray(rgb2gray(img));
df = tb(strcmp(tb.sample_id,fname) & strcmp(tb.section, 'section_2'), :);
x = df.pxl_col_in_fullres; y = df.pxl_row_in_fullres;
k = convhull(y,x);
Bmask = poly2mask(x(k), y(k), size(BW, 1), size(BW, 2)); % Black out everything outside the polygon
BmaskImg = bsxfun(@times, BW, cast(Bmask, class(BW)));
img1(~BmaskImg) = 0;

BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;
nm = NM_mask(NM_mask>0);

%Pi = sum(nm);
%maxNM = max(nm);
%Mt = size(nm,1)*maxNM;
%NMscore2 = Pi/Mt;

%Mdata.NMscore2(Mdata.sample_id == fname) = NMscore2;

HE = mean(BG_mask(:));
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
BGNMscore2 = Pi/Mt;
temp2 = sum(nm-HE < 0);

%Mdata.BGNMscore2(Mdata.sample_id == fname) = BGNMscore2;
Mdata.tempsize2(Mdata.sample_id == fname) = temp2;
end

cluster = fullfile(pwd, '/processed-data/LC_spotAndPixel_coords_25hdg75svg_louv1.txt');
clus = readtable(cluster);
clus.section = cellfun(@(x) x(end-1:end), clus.sample_id, 'UniformOutput', false);
clus.sample_id = cellfun(@(x) x(1:end-2), clus.sample_id, 'UniformOutput', false);
naIdx = strcmp(clus.section, 'NA');
naIndices = find(naIdx);
if size(naIndices, 1)> 0 , disp([naIndices, 'donot have clusters']); end

for i= 39:numel(myfiles)
fname = myfiles(i).name(1:end-4);
disp(fname);
img = imread([pwd, '/raw-data/Images/',fname,'.tif']);

NMseg_dir = fullfile(pwd, '/processed-data/Images/NMseg/');
load([NMseg_dir, fname, 'NMseg_clean.mat'])

img1 = imcomplement(mat2gray(rgb2gray(img)));
BW = img1 < 0.1;     

df = clus(strcmp(clus.sample_id,fname) & strcmp(clus.section, 's1'), :);
x = df.pxl_col_in_fullres; y = df.pxl_row_in_fullres;
k = convhull(y,x);
Bmask = poly2mask(x(k), y(k), size(BW, 1), size(BW, 2)); % Black out everything outside the polygon
BmaskImg = bsxfun(@times, BW, cast(Bmask, class(BW)));
img1(~BmaskImg) = 0;

BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;
nm = NM_mask(NM_mask>0);
HE = mean(BG_mask(:));
temp1 = sum(nm-HE < 0);
Mdata.LCtempsize1(Mdata.sample_id == fname) = temp1;

img1 = imcomplement(mat2gray(rgb2gray(img)));
BW = img1 < 0.1;     

df = clus(strcmp(clus.sample_id,fname) & strcmp(clus.section, 's2'), :);
x = df.pxl_col_in_fullres; y = df.pxl_row_in_fullres;
k = convhull(y,x);
Bmask = poly2mask(x(k), y(k), size(BW, 1), size(BW, 2)); % Black out everything outside the polygon
BmaskImg = bsxfun(@times, BW, cast(Bmask, class(BW)));
img1(~BmaskImg) = 0;

BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;
nm = NM_mask(NM_mask>0);
HE = mean(BG_mask(:));
temp2 = sum(nm-HE < 0);
Mdata.LCtempsize2(Mdata.sample_id == fname) = temp2;
end


%%
% LC cluster based BG
for i= 1:numel(myfiles)
      if ismember(i, [6, 7, 11, 18, 25, 37,38])
        continue
      end
fname = myfiles(i).name(1:end-4);
disp(fname);
img = imread([pwd, '/raw-data/Images/',fname,'.tif']);

NMseg_dir = fullfile(pwd, '/processed-data/Images/NMseg/');
load([NMseg_dir, fname, 'NMseg_clean.mat'])

img1 = imcomplement(mat2gray(rgb2gray(img)));
BW = img1 < 0.1;     

df = clus(strcmp(clus.sample_id,fname) & strcmp(clus.section, 's1'), :);
x = df.pxl_col_in_fullres; y = df.pxl_row_in_fullres;
k = convhull(y,x);
Bmask = poly2mask(x(k), y(k), size(BW, 1), size(BW, 2)); % Black out everything outside the polygon
%BmaskImg = bsxfun(@times, BW, cast(Bmask, class(BW)));
img1(~Bmask) = 0;
BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;

HE = mean(BG_mask(:));
nm = NM_mask(NM_mask>0);
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
LCBGNMscore1 = Pi/Mt;

Mdata.LCBGNMscore1(Mdata.sample_id == fname) = LCBGNMscore1;

img1 = imcomplement(mat2gray(rgb2gray(img)));
BW = img1 < 0.1;     

df = clus(strcmp(clus.sample_id,fname) & strcmp(clus.section, 's2'), :);
x = df.pxl_col_in_fullres; y = df.pxl_row_in_fullres;
k = convhull(y,x);
Bmask = poly2mask(x(k), y(k), size(BW, 1), size(BW, 2)); % Black out everything outside the polygon
%BmaskImg = bsxfun(@times, BW, cast(Bmask, class(BW)));
img1(~Bmask) = 0;
BG_mask = img1; NM_mask = img1;
NM_mask(~NM) = 0; BG_mask(NM) = 0;

HE = mean(BG_mask(:));
nm = NM_mask(NM_mask>0);
Pi = sum(nm-HE);
maxNM = max(nm-HE);
Mt = size(nm,1)*maxNM;
LCBGNMscore2 = Pi/Mt;

Mdata.LCBGNMscore2(Mdata.sample_id == fname) = LCBGNMscore2;

disp(i)
end

writetable(Mdata, fullfile(pwd,'processed-data/Images/NMscore_LCBG.csv'));
