function [BWclean, masks, info] = remove_fold_and_tear(BW1, Igray, umPerPx, varargin)
% Safer fold/tear removal for NM masks.
% Name-value params let you relax/tighten behavior without editing code.

% ----- defaults (conservative) -----
params = inputParser;
addParameter(params,'TissueSensitivity', 0.60);    % 0.45→0.60 (softer mask)
addParameter(params,'EdgeBandUm',        20);      % µm from edge considered "edge"
addParameter(params,'InteriorUm',        6);       % require ≥ this distance from edge
addParameter(params,'FoldLineLenUm',     16);      % line length for fold detection
addParameter(params,'FoldMinAreaUm2',    3000);    % area gate for folds
addParameter(params,'FoldMinAR',         4.0);     % aspect ratio gate
addParameter(params,'FoldMaxSolidity',   0.90);    % concavity gate
addParameter(params,'TearMinAreaUm2',    12000);   % tears must be big
addParameter(params,'TinyLineUm2',       10^2);    % remove tiny line hits (area)
parse(params,varargin{:});
P = params.Results;

% ----- unit helpers -----
pxPerUm = 1/umPerPx;
um2px  = @(um)  max(1, round(um  * pxPerUm));
um2px2 = @(a)   max(1, round(a   * pxPerUm^2));

I   = mat2gray(Igray);
BW1 = logical(BW1);

% ---------- 1) Tissue mask (softer & more tolerant) ----------
Ig  = imgaussfilt(I, 1.0);                                   % mild smooth
T   = imbinarize(Ig,'adaptive','Sensitivity',P.TissueSensitivity);
T   = imfill(T,'holes');
T   = bwareaopen(T, um2px2(80*80));                           % drop tiny islands

% ---------- 2) Edge distances ----------
D         = bwdist(~T);
edgeBand  = D <  um2px(P.EdgeBandUm);
Interior  = D >= um2px(P.InteriorUm);

% ---------- 3) Fold candidates (line-like near edge + shape gate) ----------
L = false(size(BW1));
for ang = -15:15:15
    L = L | imopen(BW1, strel('line', um2px(P.FoldLineLenUm), ang));
end
L = logical(L) & edgeBand;
L = bwareaopen(L, um2px2(P.TinyLineUm2));                      % tiny line hits out

CC   = bwconncomp(L);
Fold = false(size(BW1));
if CC.NumObjects > 0
    R  = regionprops(CC,'Area','MajorAxisLength','MinorAxisLength','Solidity');
    AR = arrayfun(@(r) r.MajorAxisLength / max(r.MinorAxisLength,1), R);
    isFold = ([R.Area] >= um2px2(P.FoldMinAreaUm2)) & ...
             (AR >= P.FoldMinAR) & ...
             ([R.Solidity] <= P.FoldMaxSolidity);
    if any(isFold)
        Fold(vertcat(CC.PixelIdxList{isFold})) = true;
    end
end

% ---------- 4) Tears (large interior holes) ----------
Tfilled = imfill(T,'holes');
Tears   = Tfilled & ~T;
Tears   = bwareaopen(Tears, um2px2(P.TearMinAreaUm2));

% ---------- 5) Combine ----------
BW2 = BW1 & T & Interior & ~Fold & ~Tears;

% ----- guard-rails: if we nuked too much, relax in order -----
area0 = nnz(BW1); area2 = nnz(BW2);
relaxed = struct('Interior',false,'Fold',false,'EdgeBand',false);

if area0>0 && area2 < 0.25*area0
    % (a) drop the interior constraint first
    BW2 = BW1 & T & ~Fold & ~Tears;
    relaxed.Interior = true;
end
if nnz(BW2) < 0.25*area0
    % (b) relax fold removal
    BW2 = BW1 & T & ~Tears;
    relaxed.Fold = true;
end
if nnz(BW2) < 0.25*area0
    % (c) shrink edge band effect by half (we already removed Fold; use Interior lightly)
    Interior = D >= max(1, round(um2px(P.InteriorUm)/2));
    BW2 = BW1 & T & Interior & ~Tears;
    relaxed.EdgeBand = true;
end
if nnz(BW2) < 0.05*area0
    % (d) ultimate fallback: keep only T
    BW2 = BW1 & T;
end

% ---------- 6) Final polish: opening-by-reconstruction ----------
seRec  = strel('disk', um2px(3));
E      = imerode(BW2, seRec);
BWclean = imreconstruct(E, BW2);

% ---------- debug info ----------
info = struct;
info.areas = struct('BW1',area0,'T',nnz(T),'Interior',nnz(Interior), ...
                    'AfterCombine',nnz(BW2),'Final',nnz(BWclean));
info.relaxed = relaxed;

if nargout > 1
    masks = struct('T',T,'D',D,'edgeBand',edgeBand,'Interior',Interior, ...
                   'Fold',Fold,'Tears',Tears,'LineHits',L);
end
end
