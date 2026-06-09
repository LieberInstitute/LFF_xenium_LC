Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC';
	id = '/processed-data/xenium_imageProcessing_new/registrations';
	myfiles = dir(fullfile(Md,id,'*.xml'));

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
		
		HE = imread(fullfile(Md,id,[brain,'_registered.png']));
		BW = imread(fullfile(Md,id,brain,'NMnofolds.png'));
		BW_reg = imwarp(BW, T, 'OutputView', Rout);
		temp = regionprops(BW_reg, mat2gray(rgb2gray(HE)), "area", "centroid", "boundingbox", "MeanIntensity", "MaxIntensity");			
			Tprops = struct2table(temp);

			% Optional: split centroid and bbox into separate columns
			if ~isempty(Tprops)
			    Tprops.Centroid_X = Tprops.Centroid(:,1);
			    Tprops.Centroid_Y = Tprops.Centroid(:,2);
			    Tprops.BBox_X      = Tprops.BoundingBox(:,1);
			    Tprops.BBox_Y      = Tprops.BoundingBox(:,2);
			    Tprops.BBox_Width  = Tprops.BoundingBox(:,3);
			    Tprops.BBox_Height = Tprops.BoundingBox(:,4);
			    Tprops.Centroid = [];
			    Tprops.BoundingBox = [];
			end

			outcsv = fullfile(Md, id, brain, 'NM_regionprops.csv');
			writetable(Tprops, outcsv);
				
		disp(['done ',brain])	
			
	end
