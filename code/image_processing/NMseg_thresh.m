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
	BW(Ie>0.05 & Ie<0.6) = 1;
	
	cleaned = imopen(BW, strel('disk', 3)); 
	cleanedN = imclose(cleaned, strel('disk', 3)); 
	BW1 = bwareaopen(cleanedN, 500);
	stats = regionprops(logical(BW1), 'Area'); 
	
	save(fullfile(Md,brnum,'NMthresh.mat'), "stats", "BW1")
		
	perim = bwperim(BW1);
	out = uint8(imoverlay(NM, perim, [0 1 0]));
	imwrite(imresize(out, 0.5), fullfile(Md,brnum,'Nmthresh.png'))
			
	
end