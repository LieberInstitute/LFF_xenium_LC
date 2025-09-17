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

# ---- MICRO-REFINEMENT (optional but helps the last 1–2 px) ----
# Works on a cropped ROI so it's fast. Tries small subpixel tweaks to (dx,dy),
# and optionally tiny scale/angle. It maximizes Dice on the ROI.


# 1) define a ROI around tissue to avoid edges
y0r, y1r, x0r, xr1 = bbox((dapi_can>0.5), pad=256)
he_roi   = he_nuc_can[y0r:y1r, x0r:xr1]
dapi_roi = (dapi_can>0.5)[y0r:y1r, x0r:xr1]
Hr, Wr   = he_roi.shape
cx, cy   = Wr/2.0, Hr/2.0

def dice_score_after(dx_, dy_, dtheta_deg=0.0, scale_=1.0):
    # build tiny similarity transform around the ROI center
    t = AffineTransform(translation=(-cx, -cy))
    t += AffineTransform(scale=(scale_, scale_))
    t += AffineTransform(rotation=np.deg2rad(dtheta_deg))
    t += AffineTransform(translation=(cx + dx_, cy + dy_))
    moved = warp(he_roi, t.inverse, order=0, mode='constant', cval=0.0, preserve_range=True) > 0.5
    return dice(moved, dapi_roi)

# 2) search small neighborhoods; keep it light
best = {'dx': dx, 'dy': dy, 'th': 0.0, 'sc': 1.0, 'score': -1.0}

# (a) subpixel translation refine (no scale/angle)
for DY in np.linspace(dy-1.5, dy+1.5, 13):     # step 0.25 px
    for DX in np.linspace(dx-1.5, dx+1.5, 13):
        sc = dice_score_after(DX - dx + 0, DY - dy + 0, 0.0, 1.0)
        if sc > best['score']:
            best.update({'dx': DX, 'dy': DY, 'th': 0.0, 'sc': 1.0, 'score': sc})

# (b) very small angle/scale polish (optional; comment out if you don't want it)
for TH in np.linspace(-0.15, 0.15, 7):         # ±0.15°
    for SC in np.linspace(0.997, 1.003, 7):    # ±0.3% isotropic
        sc = dice_score_after(best['dx']-dx, best['dy']-dy, TH, SC)
        if sc > best['score']:
            best.update({'th': TH, 'sc': SC, 'score': sc})

# 3) report and apply to the FULL canvas
print(f"Micro-refine => dy={best['dy']:.3f}, dx={best['dx']:.3f}, "
      f"dθ={best['th']:.3f}°, s={best['sc']:.5f} (Dice_ROI={best['score']:.4f})")

# Build final transform for full-size masks/images
# Start from your existing rotation+centering result (he_nuc_can, he_can),
# then apply similarity about the FULL canvas center.
Cxf, Cyf = W/2.0, H/2.0
t_full = AffineTransform(translation=(-Cxf, -Cyf))
t_full += AffineTransform(scale=(best['sc'], best['sc']))
t_full += AffineTransform(rotation=np.deg2rad(best['th']))
t_full += AffineTransform(translation=(Cxf + best['dx'], Cyf + best['dy']))

# Recompute outputs using this single transform
henuc_reg = warp(he_nuc_can, t_full.inverse, order=0, mode='constant', cval=0.0,
                 preserve_range=True) > 0.5
henuc_reg = henuc_reg.astype(np.uint8)

he_reg = np.empty_like(he_can, dtype=np.float32)
for ch in range(he_can.shape[2]):
    he_reg[..., ch] = warp(he_can[..., ch], t_full.inverse, order=1, mode='constant',
                           cval=0.0, preserve_range=True).astype(np.float32)

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
