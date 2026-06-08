function Nmseg_removeFolds(brnum, K)

disp(brnum)
Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing_new/registrations'
fname = fullfile(Md,brnum,'Kmeans_cluster.mat')	
	
	load(fullfile(Md,brnum,'NMthresh.mat'),  "BW1")
	S = load(fname);  
	C = S.cluster;                       
    NM = C{K};                            
    Ie = rgb2gray(NM);
    Ie = mat2gray(Ie);   
	Ie(Ie<0.1) = 0;
	Ie(Ie>0.5) = 0; 
	
	initial_mask = imbinarize(Ie, 'adaptive');
	cleaned_mask = bwareaopen(initial_mask, 100);
	folds_mask = bwareaopen(cleaned_mask, 30000);
	BWrg = region_grow_intensity(Ie, folds_mask, 0.01, 0.5);
	BW_nofolds = BW1;
	BW_nofolds(BWrg==1) = 0;
	save(fullfile(Md,brnum,'NMnofolds.mat'), "BW_nofolds")
	
	imwrite(BW_nofolds, fullfile(Md,brnum,'NMnofolds.png'))	
	
end

