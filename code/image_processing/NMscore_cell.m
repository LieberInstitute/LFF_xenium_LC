function NMscore_cell(brnum)
	
	Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC';
	od = '/processed-data/xenium_imageProcessing/';

	cells = readtable(fullfile(Md, od, brnum, ['xeniumranger_NM_DAPI_', brnum], 'outs', 'cells.csv'));
	b = readtable(fullfile(Md, od, brnum, ['xeniumranger_NM_DAPI_', brnum], 'outs', 'cell_boundaries.csv')); % cell_id, x, y (vertex rows)
	img1 = imread(fullfile(Md, od, brnum, 'HE_registered_to_DAPI.tif'));  % single channel preferred
	img = mat2gray(rgb2gray(img1)); [H, W] = size(img);

    % Convert cell_id to string for reliable matching (tables often store as cellstr)
    cells_id = string(cells.cell_id);
    b_id     = string(b.cell_id);

    meanIntensity = nan(height(cells), 1);

    % Only consider boundary cell_ids that exist in cells.csv
    [isInCells, locCells] = ismember(b_id, cells_id);
    valid_ids = unique(b_id(isInCells));

    for i = 1:numel(valid_ids)
        cid = valid_ids(i);

        verts = b(b_id == cid, :);

        x = double(verts.vertex_x);
        y = double(verts.vertex_y);

        if numel(x) < 3
            continue
        end

        % clip to image bounds
        x = min(max(x, 1), W);
        y = min(max(y, 1), H);

        mask = poly2mask(x, y, H, W);

        if any(mask(:))
            m = mean(img(mask));
        else
            m = NaN;
        end

        % row in cells.csv for this cell_id
        row_idx = find(cells_id == cid, 1);
        meanIntensity(row_idx) = m;
		
		if mod(i, 500) == 0 
		    fprintf('Computed mean intensity for %d cells\n', i);
		end
    end

    % Add column and overwrite the SAME cells.csv
    cells.mean_intensity = meanIntensity;
    writetable(cells, fullfile(Md, od, brnum, ['xeniumranger_NM_DAPI_', brnum], 'outs', 'cells.csv'));

    fprintf('Wrote mean_intensity_HE into %s\n', cells_file);
end