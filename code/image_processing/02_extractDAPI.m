% Loop over all xenium-instrument folders, extract masks, save binary TIFFs
Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC';
od = '/processed-data/xenium_imageProcessing/';
in_dir = fullfile(Md, 'raw-data', 'xenium', 'xenium-instrument');

% Python imports (once)
zarr = py.importlib.import_module('zarr');
np   = py.importlib.import_module('numpy');

d = dir(in_dir);
d = d([d.isdir]);              % keep only directories
names = {d.name};
names = names(~ismember(names, {'.','..'}));

for i = 1:numel(names)
    name = names{i};
    if startsWith(name, 'old_corrupt'), continue; end

    % Extract Br number (works for BrXXXX-L, BrXXXX_re-dis, etc.)
    tok = regexp(name, '(Br\d+)', 'tokens', 'once');
    if isempty(tok), continue; end
    br = tok{1};

    % cells.zarr(.zip) path
    zarr_zip = fullfile(in_dir, name, 'cells.zarr.zip');
    zarr_dir = fullfile(in_dir, name, 'cells.zarr');
    if exist(zarr_zip, 'file')
        zpath = zarr_zip;
        store = zarr.ZipStore(zpath, pyargs('mode','r'));
    elseif exist(zarr_dir, 'dir')
        zpath = zarr_dir;
        store = zarr.DirectoryStore(zpath);
    else
        fprintf('[skip] No cells.zarr(.zip) in %s\n', name);
        continue;
    end

    try
        % Open group
        root = zarr.group(pyargs('store', store));

        % === EXACTLY like your MATLAB lines ===
        cellmask_py = root.get('masks').get(int32(1));
        nucmask_py  = root.get('masks').get(int32(0));

        % Convert to MATLAB arrays
        cellmask = double(np.array(cellmask_py));
        nucmask  = double(np.array(nucmask_py));
       
        % Output dir per brain
        out_dir = fullfile(Md, od, br);
        if ~exist(out_dir, 'dir'), mkdir(out_dir); end

        % Save TIFFs (LZW compression; widely supported)
        imwrite(cellmask,  fullfile(out_dir, 'nucmask_binary.tif'));
        imwrite(nucmask, fullfile(out_dir, 'cellmask_binary.tif'));

        fprintf('[ok] %s -> %s\n', name, br);
    catch ME
        fprintf('[error] %s: %s\n', name, ME.message);
    end

    % Close the store (important for .zip)
    try
        store.close();
    catch
    end
end

fprintf('[done]\n');
