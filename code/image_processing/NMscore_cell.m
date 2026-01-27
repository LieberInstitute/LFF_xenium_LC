cells = readtable("cells.csv");
b = readtable("cell_boundaries.csv"); % cell_id, x, y (vertex rows)
img1 = imread("intensity_image.tif");  % single channel preferred

img = mat2gray(rgb2gray(img));

cellIDs = unique(b.cell_id);
meanIntensity = nan(height(cells),1);

% Create a map from cell_id to row index in cells.csv
[isInCells, loc] = ismember(cellIDs, cells.cell_id);

for i = 1:numel(cellIDs)
    cid = cellIDs(i);
    verts = b(b.cell_id == cid, :);
    x = verts.x;
    y = verts.y;

    if numel(x) < 3
        continue
    end

    mask = poly2mask(x, y, H, W);
    if any(mask(:))
        m = mean(img(mask));
    else
        m = NaN;
    end

    if isInCells(i)
        meanIntensity(loc(i)) = m;
    end
end

cells.mean_intensity = meanIntensity;
writetable(cells, "cells_with_mean_intensity.csv");
