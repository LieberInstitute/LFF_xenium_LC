function transformBF(slide)
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
		    
		    if contains(titleAttr, 'HE')
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
		
		HE = imread(fullfile(Md,id,[brain,'_HE.png']));
		%BW = imread(fullfile(Md,id,[brain,'_nucmaskPerim.png']));
		img_registered = imwarp(HE, T, 'OutputView', Rout);
		%out=imoverlay(img_registered,BW, [0,1,0]);
		%imwrite(imresize(out,0.5), fullfile(fullfile(Md,id,[brain,'_regis_overlay.png'])))
		imwrite(img_registered, fullfile(fullfile(Md,id,[brain,'_registered.png'])))
		disp(['done ',brain])	
	end
end
%load(fullfile(Md,id,slide,[slide,'_B1.mat']), 'TMEM119');
%BW = imbinarize(TMEM119,0.02);
%perim = bwperim(BW);
%img  = imread(fullfile(Md,id,slide,[slide, '_BF_B1.tif']));
%img_registered = imwarp(img, T, 'OutputView', Rout);
%imwrite(img_registered, fullfile(Md, od, [slide, '_B1_BF_aligned.png']))
%out=imoverlay(img_registered,perim, [0,1,0]);
%imwrite(out, fullfile(Md, od, [slide, '_B1_BF_overlay.png']))
%
%load(fullfile(Md,id,slide,[slide,'_C1.mat']), 'TMEM119');
%BW = imbinarize(TMEM119,0.02);
%perim = bwperim(BW);
%img  = imread(fullfile(Md,id,slide,[slide, '_BF_C1.tif']));
%img_registered = imwarp(img, T, 'OutputView', Rout);
%imwrite(img_registered, fullfile(Md, od, [slide, '_C1_BF_aligned.png']))
%out=imoverlay(img_registered,perim, [0,1,0]);
%imwrite(out, fullfile(Md, od, [slide, '_C1_BF_overlay.png']))
%	
%load(fullfile(Md,id,slide,[slide,'_D1.mat']), 'TMEM119');
%BW = imbinarize(TMEM119,0.02);
%perim = bwperim(BW);
%img  = imread(fullfile(Md,id,slide,[slide, '_BF_D1.tif']));
%img_registered = imwarp(img, T, 'OutputView', Rout);
%imwrite(img_registered, fullfile(Md, od, [slide, '_D1_BF_aligned.png']))
%out=imoverlay(img_registered,perim, [0,1,0]);
%imwrite(out, fullfile(Md, od, [slide, '_D1_BF_overlay.png']))
