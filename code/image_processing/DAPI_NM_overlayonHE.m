Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC';
od = 'processed-data/xenium_imageProcessing';

brains = dir(fullfile(Md, od, 'Br*'));
brains = brains([brains.isdir]);

microns_per_pixel = 0.25;
scale = 0.21 / 0.25;
for i = 1:length(brains)
    brain = brains(i).name;
    fprintf('Processing %s\n', brain);

    try
        img_file = fullfile(Md, od, brain, 'HE_registered_to_DAPI.tif');
        mat_file = fullfile(Md, od, brain, 'NM_nuc_mergedLabels.mat');
        csv_file = fullfile(Md, od, brain, 'transcripts_subset.csv');
        out_file_base = fullfile(Md, od, brain, 'NmDAPI_overlay.png');
        out_dir = fullfile(Md, od, brain, 'gene_overlays');

        if ~isfile(img_file)
            fprintf('  Missing image for %s\n', brain);
            continue
        end
        if ~isfile(mat_file)
            fprintf('  Missing mergedLabels mat for %s\n', brain);
            continue
        end
        if ~isfile(csv_file)
            fprintf('  Missing transcripts_subset.csv for %s\n', brain);
            continue
        end

        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end

        % Load image and labels
        img_registered = imread(img_file);
        img = im2double(img_registered);
        temp = load(mat_file);

        if ~isfield(temp, 'mergedLabels')
            fprintf('  mergedLabels not found in %s\n', brain);
            continue
        end

        % Make nuclei perimeter overlay
        perim = bwperim(temp.mergedLabels);
        out = uint8(imoverlay(img, perim, [0 0 0]));

        % Save base overlay
        imwrite(imresize(out, 0.5), out_file_base);

       % % Read transcript table
       % T = readtable(csv_file);
       %
       % req_cols = {'x_location','y_location','feature_name'};
       % if ~all(ismember(req_cols, T.Properties.VariableNames))
       %     fprintf('  Missing required columns in %s\n', brain);
       %     continue
       % end
       %
       % % Convert to pixel coordinates
	   %	T.x_pixel = round((T.x_location / microns_per_pixel) * scale);
	   %	T.y_pixel = round((T.y_location / microns_per_pixel) * scale);
	   %	
       % % Keep coordinates inside image bounds
       % [img_h, img_w, ~] = size(out);
       % keep = T.x_pixel >= 1 & T.x_pixel <= img_w & ...
       %        T.y_pixel >= 1 & T.y_pixel <= img_h;
       % T = T(keep, :);
       %
       % % Unique genes
       % genes = unique(T.feature_name);
       %
       % for g = 1:length(genes)
       %     gene = genes{g};
       %     fprintf('    Saving overlay for %s\n', gene);
       %
       %     idx = strcmp(T.feature_name, gene);
       %
       %     if sum(idx) == 0
       %         continue
       %     end
       %
       %     f = figure('Visible', 'off');
       %     imshow(out);
       %     hold on;
       %     scatter(T.x_pixel(idx), T.y_pixel(idx), 3, 'r', 'filled');
       %     title(sprintf('%s Transcripts Overlay - %s', gene, brain), ...
       %         'Interpreter', 'none');
       %     hold off;
       %
       %     gene_safe = regexprep(gene, '[^a-zA-Z0-9_-]', '_');
       %     out_file = fullfile(out_dir, [gene_safe '_transcripts_overlay.png']);
       %     exportgraphics(f, out_file, 'Resolution', 600);
       %     close(f);
       % end

    catch ME
        fprintf('  Error processing %s: %s\n', brain, ME.message);
    end
end