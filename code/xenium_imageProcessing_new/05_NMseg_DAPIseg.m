Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC';
id = '/processed-data/xenium_imageProcessing_new';
myfiles = dir(fullfile(Md,id,'registrations','*.xml'));

for ii = 1:length(myfiles)
	xmlFile = fullfile(myfiles(ii).folder, myfiles(ii).name);
	brain = myfiles(ii).name(1:end-4);
	xDoc = xmlread(xmlFile);
	
	% Get all t2_patch nodes
	patchNodes = xDoc.getElementsByTagName('t2_patch');
	
	% Loop through patches to find the HE image
	for i = 0:patchNodes.getLength-1
	    patch = patchNodes.item(i);
	    titleAttr = char(patch.getAttribute('title'));
	    
	    if contains(titleAttr, 'sample_05')
	        % Get transform string and parse values
	        transformStr = char(patch.getAttribute('transform'));
	        tokens = regexp(transformStr, 'matrix\((.*)\)', 'tokens');
	        values = sscanf(tokens{1}{1}, '%f,');
	        a = values(1); b = values(2); c = values(3);
	        d = values(4); e = values(5); f = values(6);
		
		elseif contains(titleAttr, 'nucmaskPerim')
	        % Get target size
	        target_width = str2double(patch.getAttribute('o_width'));
	        target_height = str2double(patch.getAttribute('o_height'));
	    end
	end
	
	T = affine2d([a, b, 0; c, d, 0; e, f, 1]);
	Rout = imref2d([target_height, target_width]);
	
	%HE = imread(fullfile(Md,id,'registrations',[brain,'_HE.png']));
	load(fullfile(Md,id,brain,'DAPImasks.mat'), 'nucmask')
	BWnm = imread(fullfile(Md,id,'registrations',brain,'NMnofolds.png'));
	BWnm = BWnm > 0;
	BWnm_registered = imwarp(BWnm, T, 'OutputView', Rout) > 0;
	%BWnm_registered = imwarp(BWnm, T, 'OutputView', Rout);
	nucBW = nucmask > 0;
	BWnm_noNuc = BWnm_registered & ~nucBW;
	%HE_registered = imwarp(HE, T, 'OutputView', Rout);
	%out=imoverlay(HE_registered,bwperim(nucBW), [0,1,0]);
	%out=imoverlay(out,bwperim(BWnm_noNuc), [1,1,1]);
	%imwrite(imresize(out,0.5), fullfile(fullfile(Md,id,brain,'NmDAPI_overlay.png')))
	
	NM_labeled = uint32(bwlabel(BWnm_noNuc));
	maxNuc = uint32(max(nucmask(:)));
	idx = NM_labeled > 0;

	NM_labeled(idx) = NM_labeled(idx) + maxNuc;

	combined_mask_uint32 = uint32(nucmask);
	combined_mask_uint32(idx) = NM_labeled(idx);

	output_file = fullfile(Md, id, brain, 'combined_nucmask.npy');
	py.numpy.save(output_file, py.numpy.array(combined_mask_uint32));
	
	disp(['done ',brain])	
end
		
			
			