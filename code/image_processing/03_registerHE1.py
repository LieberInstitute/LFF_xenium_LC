#!/usr/bin/env python3
import sys, os
import numpy as np
from tifffile import imread, imwrite
from skimage import img_as_float32, exposure
from skimage.color import rgb2hed
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects, remove_small_holes, binary_opening, disk
from skimage.feature import canny
from skimage.transform import rotate, resize, AffineTransform, warp
from skimage.registration import phase_cross_correlation
from scipy.ndimage import distance_transform_edt as dist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- helpers ----------------
def dice(a, b):
    a = a.astype(bool); b = b.astype(bool)
    inter = np.count_nonzero(a & b)
    return (2.0 * inter) / (np.count_nonzero(a) + np.count_nonzero(b) + 1e-6)

def translate_slice(img, dy_i, dx_i, fill=0):
    """Fast integer translation by slicing. dy>0 down, dx>0 right."""
    H, W = img.shape[:2]
    out = np.full_like(img, fill, dtype=img.dtype)
    y0s, x0s = max(0, -dy_i), max(0, -dx_i)
    y0d, x0d = max(0,  dy_i), max(0,  dx_i)
    h = min(H - y0s, H - y0d); w = min(W - x0s, W - x0d)
    if h>0 and w>0:
        out[y0d:y0d+h, x0d:x0d+w, ...] = img[y0s:y0s+h, x0s:x0s+w, ...]
    return out

def signed_dt(mask_bool):
    # + inside tissue, - outside
    return dist(mask_bool) - dist(~mask_bool)

def bbox(mask, pad=64):
    ys, xs = np.nonzero(mask)
    if ys.size == 0:
        return (0, mask.shape[0], 0, mask.shape[1])  # fallback to full image
    y0 = max(0, ys.min() - pad); y1 = min(mask.shape[0], ys.max() + pad + 1)
    x0 = max(0, xs.min() - pad); x1 = min(mask.shape[1], xs.max() + pad + 1)
    return (y0, y1, x0, x1)

def apply_shift_skimage(img, dx, dy, order=1, cval=0.0):
    # AffineTransform uses (x, y) = (cols, rows)
    tform = AffineTransform(translation=(dx, dy))
    return warp(img, tform.inverse, order=order, mode='constant',
                cval=cval, preserve_range=True)

# ---------------- args & paths ----------------
brnum_arg  = sys.argv[1]
slide_arg  = sys.argv[2]
sample_arg = sys.argv[3]

HEPath    = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/split_samples/{slide_arg}/{sample_arg}'
DAPIPath  = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/{brnum_arg}/nucmask_binary.tif'
out_he    = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/{brnum_arg}/HE_registered_to_DAPI.tif'
out_henuc = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/{brnum_arg}/HE_nuclei_registered_to_DAPI.tif'
out_overlay  = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/{brnum_arg}/overlay_on_dapi.png'
out_overlay1 = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/{brnum_arg}/overlay_on_dapi1.png'

# ---------------- load ----------------
he   = imread(HEPath)   # RGB H&E
dapi = imread(DAPIPath) # 2D nuclei mask (0/255 or 0/1)
assert dapi.ndim == 2, f"DAPI must be 2D, got {dapi.shape}"
H, W = dapi.shape
print("HE:", he.shape, he.dtype, " DAPI:", dapi.shape, dapi.dtype)

# ---------------- H&E nuclei segmentation ----------------
hed = rgb2hed(img_as_float32(he))     # Hematoxylin/Eosin/DAB in OD space
Hc  = hed[..., 0]
Hn  = (Hc - Hc.min()) / (Hc.ptp() + 1e-8)
Hn  = exposure.equalize_adapthist(Hn, clip_limit=0.02)
t   = threshold_otsu(Hn)
he_nuc = Hn > t
he_nuc = binary_opening(he_nuc, disk(1))
he_nuc = remove_small_holes(he_nuc, area_threshold=64)
he_nuc = remove_small_objects(he_nuc, min_size=64)
he_nuc = he_nuc.astype(np.float32)

dapi_mask = (dapi > 0).astype(np.float32)

# ---------------- coarse angle search (89–91°, step 0.1) ----------------
he_edges   = canny(he_nuc, sigma=1.0).astype(np.float32)
dapi_edges = canny(dapi_mask, sigma=1.0).astype(np.float32)

target_max = 3000
ds_full = max(1, int(max(H, W) / target_max))
Hds, Wds = max(1, H // ds_full), max(1, W // ds_full)
de_ds  = resize(dapi_edges, (Hds, Wds), preserve_range=True, anti_aliasing=True).astype(np.float32)

best_angle, best_score = None, -1e9
for ang in np.arange(89.5, 90.5, 0.1):
    her = rotate(he_edges, angle=ang, resize=True, order=1, mode='constant', cval=0, preserve_range=True)
    Hr, Wr = her.shape
    ys = max(0, (Hr - H) // 2); xs = max(0, (Wr - W) // 2)
    yd = max(0, (H  - Hr) // 2); xd = max(0, (W  - Wr) // 2)
    her_can = np.zeros((H, W), dtype=np.float32)
    ys0, xs0 = ys, xs
    ye0, xe0 = ys + min(H, Hr), xs + min(W, Wr)
    yd0, xd0 = yd, xd
    her_can[yd0:yd0+(ye0-ys0), xd0:xd0+(xe0-xs0)] = her[ys0:ye0, xs0:xe0]
    her_ds = resize(her_can, (Hds, Wds), preserve_range=True, anti_aliasing=True).astype(np.float32)
    a = (her_ds - her_ds.mean()) / (her_ds.std() + 1e-6)
    b = (de_ds  - de_ds.mean())  / (de_ds.std()  + 1e-6)
    score = float((a * b).mean())
    if score > best_score:
        best_score, best_angle = score, float(ang)

print(f"Best angle ≈ {best_angle:.2f}° (score={best_score:.4f})")

# ---------------- rotate & center to DAPI canvas ----------------
# nuclei (nearest)
he_nuc_rot = rotate(he_nuc, angle=best_angle, resize=True,
                    order=0, mode='constant', cval=0, preserve_range=True).astype(np.float32)
Hr, Wr = he_nuc_rot.shape
ys = max(0, (Hr - H) // 2); xs = max(0, (Wr - W) // 2)
yd = max(0, (H  - Hr) // 2); xd = max(0, (W  - Wr) // 2)
he_nuc_can = np.zeros((H, W), dtype=np.float32)
he_nuc_can[yd:yd+min(H, Hr), xd:xd+min(W, Wr)] = he_nuc_rot[ys:ys+min(H, Hr), xs:xs+min(W, Wr)]

dapi_can = dapi_mask.astype(np.float32)

# quick pre-translation QC
plt.figure(figsize=(8,8))
plt.imshow(dapi_can, cmap='gray'); plt.contour(he_nuc_can > 0.5, levels=[0.5], colors=['magenta'], linewidths=0.7)
plt.axis('off'); plt.tight_layout(); plt.savefig(out_overlay1, dpi=300, bbox_inches='tight', pad_inches=0); plt.close()
print("Saved overlay before translation:", out_overlay1)

# ---------------- robust translation estimation ----------------
he_bin   = (he_nuc_can > 0.5)
dapi_bin = (dapi_can   > 0.5)

# 1) crop to tissue bbox (avoid blank paper/canvas)
y0, y1, x0, x1 = bbox(dapi_bin)
he_cut   = he_bin[y0:y1, x0:x1]
dapi_cut = dapi_bin[y0:y1, x0:x1]
Hc, Wc = dapi_cut.shape

# 2) downsample + signed distance transforms
target_max = 3000
ds = max(1, int(max(Hc, Wc) / target_max))
Hds2, Wds2 = max(1, Hc // ds), max(1, Wc // ds)

he_ds   = resize(he_cut.astype(float),   (Hds2, Wds2), preserve_range=True) > 0.5
dapi_ds = resize(dapi_cut.astype(float), (Hds2, Wds2), preserve_range=True) > 0.5

heSDT_ds   = signed_dt(he_ds).astype(np.float32)
dapiSDT_ds = signed_dt(dapi_ds).astype(np.float32)

wy = np.hanning(Hds2); wx = np.hanning(Wds2); win = (wy[:, None] * wx[None, :]).astype(np.float32)
heSDT_ds_w   = heSDT_ds * win
dapiSDT_ds_w = dapiSDT_ds * win

# 3) initial shift (DS space), then convert to full-res
shift, _, _ = phase_cross_correlation(dapiSDT_ds_w, heSDT_ds_w, upsample_factor=100)
dy0_ds, dx0_ds = float(shift[0]), float(shift[1])           # rows, cols (DS grid)

# 4) sign check + coarse-to-fine Dice refinement (DS grid)
def dice_at_ds(dx_ds, dy_ds):
    test = translate_slice(he_ds.astype(np.uint8), int(round(dy_ds)), int(round(dx_ds)), 0)
    return dice(test, dapi_ds)

candidates = [(dx0_ds, dy0_ds), (-dx0_ds, -dy0_ds)]
best_dx_ds, best_dy_ds, best_score = 0.0, 0.0, -1.0
for DX, DY in candidates:
    sc = dice_at_ds(DX, DY)
    if sc > best_score:
        best_dx_ds, best_dy_ds, best_score = DX, DY, sc

for step, rad in [(8, 128), (2, 24)]:
    dx_grid = np.arange(best_dx_ds - rad, best_dx_ds + rad + 1e-6, step)
    dy_grid = np.arange(best_dy_ds - rad, best_dy_ds + rad + 1e-6, step)
    for DY in dy_grid:
        for DX in dx_grid:
            sc = dice_at_ds(DX, DY)
            if sc > best_score:
                best_dx_ds, best_dy_ds, best_score = DX, DY, sc

dx, dy = best_dx_ds * ds, best_dy_ds * ds  # cols, rows in FULL pixels
print(f"Final translation picked by Dice: dy={dy:.2f}, dx={dx:.2f} (Dice@DS={best_score:.3f})")

# ---------------- apply translation (skimage.warp) ----------------
# rotate HE RGB once (already done for nuclei)
he_f   = img_as_float32(he)
he_rot = rotate(he_f, angle=best_angle, resize=True,
                order=1, mode='constant', cval=0, preserve_range=True).astype(np.float32)
Hr, Wr = he_rot.shape[:2]
ys = max(0, (Hr - H) // 2); xs = max(0, (Wr - W) // 2)
yd = max(0, (H  - Hr) // 2); xd = max(0, (W  - Wr) // 2)
he_can = np.zeros((H, W, he_rot.shape[2]), dtype=np.float32)
he_can[yd:yd+min(H, Hr), xd:xd+min(W, Wr), :] = he_rot[ys:ys+min(H, Hr), xs:xs+min(W, Wr), :]

# nuclei (nearest)
henuc_reg = apply_shift_skimage(he_nuc_can, dx, dy, order=0) > 0.5
henuc_reg = henuc_reg.astype(np.uint8)

# RGB (bilinear)
he_reg = np.empty_like(he_can, dtype=np.float32)
for ch in range(he_can.shape[2]):
    he_reg[..., ch] = apply_shift_skimage(he_can[..., ch], dx, dy, order=1, cval=0.0).astype(np.float32)

# ---------------- save ----------------
imwrite(out_henuc, henuc_reg)
imwrite(out_he, (np.clip(he_reg, 0, 1) * 65535).astype(np.uint16))
print("Wrote:", out_henuc)
print("Wrote:", out_he)

# ---------------- QC overlay ----------------
plt.figure(figsize=(8,8))
plt.imshow(np.clip(he_reg, 0, 1))
plt.contour(dapi_mask > 0, levels=[0.5], colors=['lime'], linewidths=0.7)
plt.axis('off'); plt.tight_layout(); plt.savefig(out_overlay, dpi=300, bbox_inches='tight', pad_inches=0); plt.close()
print("Wrote:", out_overlay)
