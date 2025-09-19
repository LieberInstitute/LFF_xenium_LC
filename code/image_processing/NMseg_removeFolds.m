function Nmseg_removeFolds(brnum, K)

disp(brnum)
Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/'
fname = fullfile(Md,brnum,'Kmeans_cluster.mat')	
	S = load(fname);  
	C = S.cluster;                       
    NM = C{K};                            
    Ie = rgb2gray(NM);
    Ie = mat2gray(Ie);   
	Ie(Ie<0.1) = 0;
	Ie(Ie>0.5) = 0; 
	
	initial_mask = imbinarize(gray_image, 'adaptive');
	cleaned_mask = bwareaopen(initial_mask, 100);
	folds_mask = bwareaopen(cleaned_mask, 25000);
	BWrg = region_grow_intensity(Ie, folds_mask, 0.01, 0.5);
	
	
		imwrite(BWrg, fullfile(Md,brnum,'Nmthresh_nofolds.png'))	
	
end
		
		%%% testing
				
%bad = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/Br1039-L/Kmeans_cluster.mat';  
%good = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/Br0946/Kmeans_cluster.mat'; 
%load(bad)
%	B = cluster{1};
%load(good)   
%	G = cluster{1};
%
%GIe = mat2gray(rgb2gray(G));
%BIe = mat2gray(rgb2gray(B));
%
%BWG = zeros(size(GIe)); BWG(GIe>0.05 & GIe<graythresh(GIe)) =1;
%BWB = zeros(size(BIe)); BWB(BIe>0.05 & BIe<graythresh(BIe)) =1;
%
%Bc = imopen(BWB, strel('disk', 3)); 
%BcN = imclose(Bc, strel('disk', 3)); 
%bwB = bwareaopen(BcN, 500);
%
%Gc = imopen(BWG, strel('disk', 3)); 
%GcN = imclose(Gc, strel('disk', 3)); 
%bwG = bwareaopen(GcN, 500);
%

