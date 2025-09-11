Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
od ='/processed-data/xenium_imageProcessing/';
brain = 'Br6423';
imgname = [brain, '_HE_21.png'];

xmlFile = fullfile(Md, od,brain, [brain, '.xml']);
xDoc = xmlread(xmlFile);

% Get all t2_patch nodes
patchNodes = xDoc.getElementsByTagName('t2_patch');

% Loop through patches to find the HE image
for i = 0:patchNodes.getLength-1
    patch = patchNodes.item(i);
    titleAttr = char(patch.getAttribute('title'));
    
    if contains(titleAttr, imgname)
        % Get transform string and parse values
        transformStr = char(patch.getAttribute('transform'));
        tokens = regexp(transformStr, 'matrix\((.*)\)', 'tokens');
        values = sscanf(tokens{1}{1}, '%f,');
        a = values(1); b = values(2); c = values(3);
        d = values(4); e = values(5); f = values(6);
	
	elseif contains(titleAttr, 'cellmask.png')
        % Get target size
        target_width = str2double(patch.getAttribute('o_width'));
        target_height = str2double(patch.getAttribute('o_height'));
    end
end

img = imread(fullfile(Md,od,brain,imgname));
T = affine2d([a, b, 0; c, d, 0; e, f, 1]);
Rout = imref2d([target_height, target_width]);
img_registered = imwarp(img, T, 'OutputView', Rout);

imwrite(img_registered, fullfile(Md, od, brain, [brain, '_HE_aligned.png']))