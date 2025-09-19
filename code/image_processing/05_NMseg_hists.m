txt = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/code/image_processing/inputs.txt';

% Build robust options for a tab-separated file
opts = delimitedTextImportOptions('Delimiter', '\t', ...
    'DataLines', 1, 'NumVariables', 5);

opts.VariableNames = {'Br', 'Slide', 'File', 'Rank', 'Note'};
opts.VariableTypes = {'string','string','string','double','string'};
opts.ExtraColumnsRule = 'ignore';
opts.EmptyLineRule = 'read';
opts = setvaropts(opts, {'Br','Slide','File','Note'}, ...
                  'WhitespaceRule','preserve', 'EmptyFieldRule','auto');

T = readtable(txt, opts);

MD = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing';
n = height(T);
cmap = lines(n);   xgrid = linspace(0,1,256);  
edges = linspace(0,1,256);

figure; hold on
for i = 1:n
    br = T.Br(i);                                   % string
    matfile = fullfile(MD, br, 'Kmeans_cluster.mat');
    load(matfile);  
 r = T.Rank(i);                               
    NM = cluster{r};                            
    Ie = rgb2gray(NM);
    Ie = mat2gray(Ie);                         
 
 histogram(Ie(:), ...
         'BinEdges', edges, ...
         'Normalization','probability', ...
         'DisplayStyle','stairs', ...
         'EdgeColor', cmap(i,:), ...
         'LineWidth', 1.0);
			  
%[f,xi] = ksdensity(Ie(:), xgrid);
%plot(xi, f, 'Color', cmap(i,:), 'LineWidth', 1.2);
clear cluster
 disp(i)
end
hold off
xlabel('Intensity'); ylabel('Probability');
title('NM class intensity histograms across brains');
