import sys,os
os.chdir('/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/')
import warnings
import numpy as np
import tifffile as tiff
import cv2
from skimage import morphology, measure
from pyhere import here

slide_dir = here("raw-data/xenium/post-xenium_images")
fname_arg = sys.argv[1]
fname = os.path.join(slide_dir, fname_arg)

# derive output directory based on sample ID
basename = os.path.basename(fname)   # e.g. LC_PostXenium_0068968_Slide5_40x_05_13_2025_HRD.tif
parts = basename.split("_")
sample_id = "_".join(parts[2:4])     # "0068968_Slide5"
outdir = here("processed-data/xenium_imageProcessing/split_samples", sample_id)
os.makedirs(outdir, exist_ok=True)

print("Input file :", fname)
print("Output dir :", outdir)

# Optional tuning knobs:
PREVIEW_MAXDIM = 4000     # preview max dimension (pixels)
MARGIN_FULLRES = 200      # padding around each box at full-res
MIN_AREA_FRAC = 0.0005    # min component area relative to preview image area
EXPECTED_N = 7            # if >0, keep largest N components; set 0 to keep all

# ---- preview reader (pick smallest pyramid) ----
def read_preview_with_tifffile(path, maxdim=4000):
    """Return a low-res RGB preview and scaling factors."""
    with tiff.TiffFile(path) as tf:
        series = tf.series[0]
        # try to use the smallest-resolution page
        try:
            pages = series.pages
            areas = [p.shape[0]*p.shape[1] for p in pages if len(p.shape) >= 2]
            min_idx = int(np.argmin(areas))
            preview = pages[min_idx].asarray()
            full_h, full_w = series.pages[0].shape[:2]
        except Exception:
            # fallback: just read first page
            preview = series.pages[0].asarray()
            full_h, full_w = preview.shape[:2]
    # force to RGB uint8
    if preview.ndim == 2:
        preview = np.stack([preview]*3, axis=-1)
    if preview.dtype != np.uint8:
        pmin, pmax = np.percentile(preview, (1, 99))
        if pmax <= pmin:
            pmin, pmax = preview.min(), preview.max()
        preview = np.clip((preview - pmin) / max(pmax - pmin, 1e-6), 0, 1)
        preview = (preview * 255).astype(np.uint8)
    h, w = preview.shape[:2]
    scale = min(maxdim / max(w, h), 1.0)
    if scale < 1.0:
        new_w = int(round(w * scale))
        new_h = int(round(h * scale))
        preview = cv2.resize(preview, (new_w, new_h), interpolation=cv2.INTER_AREA)
        sx = full_w / new_w
        sy = full_h / new_h
    else:
        sx, sy = full_w / w, full_h / h
    return preview, full_w, full_h, sx, sy
    
# ---- tissue mask on preview ----
def tissue_mask_from_preview(rgb_u8):
    hsv = cv2.cvtColor(rgb_u8, cv2.COLOR_RGB2HSV)
    h, s, v = cv2.split(hsv)
    sat_mask = s > 20
    inv_gray = 255 - cv2.cvtColor(rgb_u8, cv2.COLOR_RGB2GRAY)
    _, otsu = cv2.threshold(inv_gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    mask = (otsu > 0) & sat_mask
    # clean up
    mask = morphology.remove_small_holes(mask, area_threshold=500)
    mask = morphology.remove_small_objects(mask, min_size=500)
    mask = morphology.binary_closing(mask, morphology.disk(5))
    mask = morphology.binary_opening(mask, morphology.disk(3))
    return mask.astype(np.uint8)
    
    
# ---- find components & boxes ----
def find_components(mask, min_area_frac=0.0005):
    labels = measure.label(mask, connectivity=2)
    props = measure.regionprops(labels)
    H, W = mask.shape
    thr = min_area_frac * (H*W)
    keep = [p for p in props if p.area >= thr]
    return keep
    
def expand_box(y0, x0, y1, x1, pad, H, W):
    return max(0, y0-pad), max(0, x0-pad), min(H, y1+pad), min(W, x1+pad)
    
# ---- region reader: lazy if possible ----
def read_region_lazy(path, y0, y1, x0, x1):
    try:
        import zarr
        store = tiff.imread(path, aszarr=True)
        za = zarr.open(store, mode="r")
        if za.ndim == 2:
            region = za[y0:y1, x0:x1]
        elif za.ndim == 3:
            if za.shape[2] in (1, 3, 4):  # HWC
                region = za[y0:y1, x0:x1, :]
            else:  # CHW -> HWC
                region = np.moveaxis(za, 0, -1)[y0:y1, x0:x1, :]
        else:
            raise RuntimeError(f"Unsupported zarr dims: {za.shape}")
        return np.array(region)
    except Exception:
        warnings.warn("Lazy zarr read unavailable; loading full image (may be large).")
        arr = tiff.imread(path)
        if arr.ndim == 2:
            return arr[y0:y1, x0:x1]
        elif arr.ndim == 3:
            if arr.shape[-1] in (1, 3, 4):  # HWC
                return arr[y0:y1, x0:x1, :]
            else:  # CHW -> HWC
                return np.moveaxis(arr, 0, -1)[y0:y1, x0:x1, :]
                
# ---- main ----
print(f"[INFO] Reading preview: {fname}")
preview, full_w, full_h, sx, sy = read_preview_with_tifffile(fname, PREVIEW_MAXDIM)
print(f"[INFO] Full-res: {full_w} x {full_h} | Preview: {preview.shape[1]} x {preview.shape[0]} | scale sx={sx:.3f}, sy={sy:.3f}")

mask = tissue_mask_from_preview(preview)
props = find_components(mask, min_area_frac=MIN_AREA_FRAC)

if EXPECTED_N and len(props) > EXPECTED_N:
    # keep the largest N by area
    props = sorted(props, key=lambda r: r.area, reverse=True)[:EXPECTED_N]

# sort left->right, then top->bottom for reproducible order
props = sorted(props, key=lambda p: (p.bbox[1], p.bbox[0]))

if not props:
    raise SystemExit("[ERROR] No tissue components found. Tune MIN_AREA_FRAC or check the slide.")

print(f"[INFO] Found {len(props)} tissue blobs. Exporting crops to: {outdir}")

for i, p in enumerate(props, 1):
    min_row, min_col, max_row, max_col = p.bbox  # preview coords
    # map to full-res
    y0 = int(round(min_row * sy))
    x0 = int(round(min_col * sx))
    y1 = int(round(max_row * sy))
    x1 = int(round(max_col * sx))
    y0, x0, y1, x1 = expand_box(y0, x0, y1, x1, MARGIN_FULLRES, full_h, full_w)
    region = read_region_lazy(fname, y0, y1, x0, x1)
    # ensure 8-bit RGB for portability when saving
    arr = region
    if arr.dtype != np.uint8:
        pmin, pmax = np.percentile(arr, (1, 99))
        if pmax <= pmin:
            pmin, pmax = float(arr.min()), float(arr.max())
        arr = np.clip((arr - pmin) / max(pmax - pmin, 1e-8), 0, 1)
        arr = (arr * 255).astype(np.uint8)
    if arr.ndim == 2:
        arr = np.stack([arr]*3, axis=-1)
    elif arr.ndim == 3 and arr.shape[2] == 1:
        arr = np.repeat(arr, 3, axis=2)
    savepath = os.path.join(outdir, f"sample_{i:02d}_y{y0}_x{x0}_h{y1-y0}_w{x1-x0}.tif")
    scale = 0.25 / 0.21
    arr_resized = cv2.resize(arr, None, fx=scale, fy=scale, interpolation=cv2.INTER_LINEAR)
    tiff.imwrite(savepath, arr_resized, photometric='rgb')
    print(f"[OK] Saved: {savepath}")
    
# ---- montage (low-res) from preview showing each split + sample number ----
import math

def make_numbered_thumb(img, label, target_h=320):
    # resize to fixed height, keep aspect
    h, w = img.shape[:2]
    scale = target_h / max(h, 1)
    tw = max(1, int(round(w * scale)))
    thumb = cv2.resize(img, (tw, target_h), interpolation=cv2.INTER_AREA)

    # draw a white box + black number in the top-left
    pad = 8
    box_w, box_h = 64, 44
    cv2.rectangle(thumb, (pad, pad), (pad + box_w, pad + box_h), (255, 255, 255), thickness=-1)
    cv2.putText(thumb, str(label), (pad + 10, pad + 32),
                cv2.FONT_HERSHEY_SIMPLEX, 1.0, (0, 0, 0), 2, cv2.LINE_AA)
    return thumb

def stack_with_gaps(row_imgs, gap=12, bg=(255, 255, 255)):
    """Horizontally stack a list of same-height images with gaps."""
    if not row_imgs:
        return None
    h = row_imgs[0].shape[0]
    # insert vertical gap blocks
    pieces = []
    for i, im in enumerate(row_imgs):
        pieces.append(im)
        if i < len(row_imgs) - 1:
            pieces.append(np.full((h, gap, 3), bg, dtype=np.uint8))
    return np.hstack(pieces)

def build_montage(thumbs, cols=4, gap=12, bg=(255, 255, 255)):
    """Arrange thumbs into a grid with gaps and return an RGB montage array."""
    if not thumbs:
        return None
    n = len(thumbs)
    cols = max(1, min(cols, n))
    rows = math.ceil(n / cols)

    # pad the last row with blank tiles to make a full grid
    target_h = thumbs[0].shape[0]
    blank = np.full((target_h, max(t.shape[1] for t in thumbs), 3), bg, dtype=np.uint8)

    # split into rows
    rows_imgs = []
    k = 0
    for r in range(rows):
        row = []
        for c in range(cols):
            if k < n:
                row.append(thumbs[k])
            else:
                row.append(blank[:, :thumbs[0].shape[1]])  # crop blank to typical width
            k += 1
        # make same height and pad widths
        h = max(im.shape[0] for im in row)
        wmax = max(im.shape[1] for im in row)
        row_padded = []
        for im in row:
            # pad right side to wmax
            pad_w = wmax - im.shape[1]
            if pad_w > 0:
                im = np.hstack([im, np.full((im.shape[0], pad_w, 3), bg, dtype=np.uint8)])
            row_padded.append(im)
        rows_imgs.append(stack_with_gaps(row_padded, gap=gap, bg=bg))

    # stack rows with horizontal gap bars
    maxw = max(rimg.shape[1] for rimg in rows_imgs)
    gap_bar = np.full((gap, maxw, 3), bg, dtype=np.uint8)
    canvas_parts = []
    for i, rimg in enumerate(rows_imgs):
        # pad row to max width
        pad_w = maxw - rimg.shape[1]
        if pad_w > 0:
            rimg = np.hstack([rimg, np.full((rimg.shape[0], pad_w, 3), bg, dtype=np.uint8)])
        canvas_parts.append(rimg)
        if i < len(rows_imgs) - 1:
            canvas_parts.append(gap_bar)
    return np.vstack(canvas_parts)

# Build thumbnails directly from preview-space crops
thumbs = []
for idx, p in enumerate(props, 1):
    r0, c0, r1, c1 = p.bbox  # preview coords
    # add a tiny pad in preview pixels so the thumb shows some margin
    pad_preview = 6
    r0 = max(0, r0 - pad_preview)
    c0 = max(0, c0 - pad_preview)
    r1 = min(preview.shape[0], r1 + pad_preview)
    c1 = min(preview.shape[1], c1 + pad_preview)
    crop_prev = preview[r0:r1, c0:c1].copy()
    thumbs.append(make_numbered_thumb(crop_prev, idx, target_h=320))

# Choose columns (e.g., up to 4 per row)
montage = build_montage(thumbs, cols=min(4, len(thumbs)), gap=12, bg=(255, 255, 255))
if montage is not None:
    montage_path = os.path.join(outdir, "montage_preview_samples.tif")  # or .png
    # Save as RGB TIFF
    tiff.imwrite(montage_path, montage, photometric='rgb')
    print(f"[OK] Saved montage: {montage_path}")


print("[DONE]")