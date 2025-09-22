function NMseg_DAPI(brnum)
	
	Md = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC';
	od = '/processed-data/xenium_imageProcessing/';
	load(fullfile(Md, od, brnum, 'DAPImasks.mat'), 'nucmask');   % expects variable: nucmask (binary or labeled)
	load(fullfile(Md, od, brnum, 'NMnofolds.mat'));   % expects variable: BW_nofolds (binary)

	%% Parameters
	minFrac = 0.10;   % 30% of nuc area
	NM_labeled = bwlabel(BW_nofolds);
	maxNuc = max(nucmask(:));
	maxNM  = max(NM_labeled(:));

	%% Precompute nucleus areas
	nucAreas = accumarray(nonzeros(nucmask), 1, [maxNuc, 1]);

	%% Output label image starts as nuclei (we keep these IDs as-is)
	mergedLabels = nucmask;

	%% We’ll assign new NM-only IDs starting after maxNuc
	nextFreeID = maxNuc + 1;

	% Optional: bookkeeping for QC
	nm_ids      = (1:maxNM).';
	best_nuc_id = zeros(maxNM,1);
	overlap_px  = zeros(maxNM,1);
	overlap_fr  = zeros(maxNM,1);
	decision    = strings(maxNM,1);  % "merged_to_nuc" or "new_nm_id"
	
	disp(maxNM)
	%% Process each NM object
	for k = 1:maxNM
	    pix = (NM_labeled == k);
	 	% Which nuclei (if any) overlap this NM object?
	    nuc_over = nucmask(pix);
	    nuc_over = nuc_over(nuc_over > 0);
	    if isempty(nuc_over)
	        % No overlap at all -> give a new ID
	        mergedLabels(pix) = nextFreeID;
	        best_nuc_id(k)  = 0;
	        overlap_px(k)   = 0;
	        overlap_fr(k)   = 0;
	        decision(k)     = "new_nm_id";
	        nextFreeID = nextFreeID + 1;
	        continue;
	    end

	    % Count overlap pixels per nucleus
	    counts = accumarray(nuc_over, 1, [maxNuc, 1]);
	    [mx, nucID] = max(counts);  % best-overlap nucleus
		
	    % Compute fraction relative to THAT nucleus’s area
			nmArea  = nnz(NM_labeled == k);
			nucArea = nucAreas(nucID);
			%frac_nm  = mx / max(nmArea, 1);   % overlap as fraction of NM
			frac_nuc = mx /nucArea;  % overlap as fraction of nucleus

			% Merge rule options:
			if frac_nuc >= minFrac 
			    mergedLabels(pix) = nucID;
			    decision(k) = "merged_to_nuc";
				best_nuc_id(k) = nucID;
				overlap_px(k) = mx;
				overlap_fr(k) = frac_nuc;
			else
			    mergedLabels(pix) = nextFreeID;
			    decision(k) = "new_nm_id";
			    nextFreeID = nextFreeID + 1;
			end
	disp(k)
	end

	%% Build a QC table
	NM_merge_summary = table(nm_ids, best_nuc_id, overlap_px, overlap_fr, decision, ...
	    'VariableNames', {'NM_ID','MatchedNucID','OverlapPixels','%NucOverlap','Decision'});

	%% (Optional) Save outputs
	outDir = fullfile(Md, od, brnum);
	if ~exist(outDir, 'dir'); mkdir(outDir); end
	save(fullfile(outDir, 'NM_nuc_mergedLabels.mat'), 'mergedLabels', 'NM_merge_summary', 'minFrac', '-v7.3');

	% If you want a quick sanity image:
	imwrite(label2rgb(mergedLabels, 'jet', 'k'), fullfile(outDir,'NM_nuc_mergedLabels_preview.png'));
	
	%% save as Numpy array
	combined_mask_uint32 = uint32(mergedLabels);
	np_array = py.numpy.array(combined_mask_uint32);

	output_file = fullfile(Md, od, brnum, 'combined_nucmask.npy');
	py.numpy.save(output_file, np_array);
	
end