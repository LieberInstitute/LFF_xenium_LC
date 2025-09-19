function NMseg_thresh(brnum, K)

disp(brnum)
Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/'
fname = fullfile(Md,brnum,'Kmeans_cluster.mat')	
	S = load(fname);  
	C = S.cluster;                       
    NM = C{K};                            
    Ie = rgb2gray(NM);
    Ie = mat2gray(Ie);   
	BW = zeros(size(Ie));
	%BW(Ie>0.05 & Ie<0.6) = 1;
	BW(Ie>0.05 & Ie<graythresh(Ie)) = 1;
	
	cleaned = imopen(BW, strel('disk', 3)); 
	cleanedN = imclose(cleaned, strel('disk', 3)); 
	BW1 = bwareaopen(cleanedN, 500);
	stats = regionprops(logical(BW1), 'Area'); 
	
	save(fullfile(Md,brnum,'NMthresh.mat'), "stats", "BW1")
		
	%perim = bwperim(BW1);
	%out = uint8(imoverlay(NM, perim, [0 1 0]));
	%imwrite(imresize(out, 0.5), fullfile(Md,brnum,'Nmthresh_overlay.png'))
	imwrite(BW1, fullfile(Md,brnum,'Nmthresh1.png'))	
	
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

