import numpy as np, cv2, matplotlib.pyplot as plt
from tifffile import imread, imwrite
from skimage import img_as_float32, exposure
from skimage.color import rgb2hed
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects, remove_small_holes, binary_opening, disk
from skimage.feature import canny
from skimage.transform import rotate, resize
from skimage.registration import phase_cross_correlation
import matplotlib
matplotlib.use("Agg")  # non-GUI backend
import matplotlib.pyplot as plt
import sys

try:
    from scipy.ndimage import distance_transform_edt as dist
except Exception:
    dist = None  # will fall back to edges

brnum_arg = sys.argv[1]
slide_arg = sys.argv[2]
sample_arg = sys.argv[3]

HEPath   = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/split_samples/{slide_arg}/{sample_arg}'
DAPIPath = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/{brnum_arg}/nucmask_binary.tif'
out_he       = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/{brnum_arg}/HE_registered_to_DAPI.tif'
out_henuc    = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/{brnum_arg}/HE_nuclei_registered_to_DAPI.tif'
out_overlay  = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/{brnum_arg}/overlay_on_dapi.png'
out_overlay1  = f'/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/{brnum_arg}/overlay_on_dapi1.png'


# ---------------- load ----------------
he   = imread(HEPath)   # RGB H&E
dapi = imread(DAPIPath) # 2D mask (cells)
assert dapi.ndim == 2, f"DAPI must be 2D, got {dapi.shape}"
H, W = dapi.shape
print("HE:", he.shape, he.dtype, " DAPI:", dapi.shape, dapi.dtype)


# ---------------- H&E nuclei segmentation ----------------
hed = rgb2hed(img_as_float32(he))     # Hematoxylin/Eosin/DAB in OD space
Hc  = hed[..., 0]
Hn  = Hc.copy()
Hn  = (Hn - np.min(Hn)) / (np.ptp(Hn) + 1e-8)
Hn  = exposure.equalize_adapthist(Hn, clip_limit=0.02)
t   = threshold_otsu(Hn)
he_nuc = Hn > t
he_nuc = binary_opening(he_nuc, disk(1))
he_nuc = remove_small_holes(he_nuc, area_threshold=64)
he_nuc = remove_small_objects(he_nuc, min_size=64)
he_nuc = he_nuc.astype(np.float32)

dapi_mask = (dapi > 0).astype(np.float32)

# ---------------- coarse angle search (89–91°, step 0.1) using edges ----------------
he_edges   = canny(he_nuc, sigma=1.0).astype(np.float32)
dapi_edges = canny(dapi_mask, sigma=1.0).astype(np.float32)


# Downsample for speed (longest side ≈ 3000 px)
target_max = 3000
ds = max(1, int(max(H, W) / target_max))
Hds, Wds = max(1, H // ds), max(1, W // ds)
de_ds  = resize(dapi_edges, (Hds, Wds), preserve_range=True, anti_aliasing=True).astype(np.float32)

best_angle, best_score = None, -1e9
for ang in np.arange(89.5, 90.5, 0.1):
    her = rotate(he_edges, angle=ang, resize=True, order=1, mode='constant', cval=0, preserve_range=True)
    # center to DAPI size (downsampled)
    Hr, Wr = her.shape
    # compute center-to-canvas (downsampled) offsets for this angle
    ys = max(0, (Hr - H) // 2); xs = max(0, (Wr - W) // 2)
    yd = max(0, (H  - Hr) // 2); xd = max(0, (W  - Wr) // 2)
    # crop/pad to (H,W)
    her_can = np.zeros((H, W), dtype=np.float32)
    ys0, xs0 = ys, xs
    ye0, xe0 = ys + min(H, Hr), xs + min(W, Wr)
    yd0, xd0 = yd, xd
    her_can[yd0:yd0+(ye0-ys0), xd0:xd0+(xe0-xs0)] = her[ys0:ye0, xs0:xe0]
    # downsample to match de_ds
    her_ds = resize(her_can, (Hds, Wds), preserve_range=True, anti_aliasing=True).astype(np.float32)
    a = (her_ds - her_ds.mean()) / (her_ds.std() + 1e-6)
    b = (de_ds  - de_ds.mean())  / (de_ds.std()  + 1e-6)
    score = float((a * b).mean())
    if score > best_score:
        best_score, best_angle = score, float(ang)

print(f"Best angle ≈ {best_angle:.2f}° (score={best_score:.4f})")

# ---------------- rotate & center both to DAPI frame (no rescaling), then translation ----------------
# Rotate he_nuc and center to DAPI canvas
he_nuc_rot = rotate(he_nuc, angle=best_angle, resize=True,
                    order=0, mode='constant', cval=0, preserve_range=True).astype(np.float32)
Hr, Wr = he_nuc_rot.shape
ys = max(0, (Hr - H) // 2); xs = max(0, (Wr - W) // 2)
yd = max(0, (H  - Hr) // 2); xd = max(0, (W  - Wr) // 2)
he_nuc_can = np.zeros((H, W), dtype=np.float32)
ys0, xs0 = ys, xs
ye0, xe0 = ys + min(H, Hr), xs + min(W, Wr)
yd0, xd0 = yd, xd
he_nuc_can[yd0:yd0+(ye0-ys0), xd0:xd0+(xe0-xs0)] = he_nuc_rot[ys0:ye0, xs0:xe0]

dapi_can = dapi_mask.astype(np.float32)

#QC
plt.figure(figsize=(8,8))
plt.imshow(dapi_mask, cmap='gray')  # DAPI background
plt.contour(he_nuc_can > 0.5, levels=[0.5], colors=['magenta'], linewidths=0.7)  # nuclei contour
plt.axis('off'); plt.tight_layout()
plt.savefig(out_overlay1, dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()

print("Saved overlay")

# Build registration fields (distance transforms preferred)
if dist is not None:
    he_field   = (dist(he_nuc_can > 0)).astype(np.float32)
    dapi_field = (dist(dapi_can   > 0)).astype(np.float32)
else:
    he_field   = canny(he_nuc_can, sigma=1.0).astype(np.float32)
    dapi_field = canny(dapi_can,   sigma=1.0).astype(np.float32)

# Hann window to reduce boundary effects
wy = np.hanning(H); wx = np.hanning(W)
win = (wy[:, None] * wx[None, :]).astype(np.float32)
heF   = he_field   * win
dapiF = dapi_field * win

# Downsample for speed (reuse ds)
heF_ds   = resize(heF,   (Hds, Wds), preserve_range=True, anti_aliasing=True).astype(np.float32)
dapiF_ds = resize(dapiF, (Hds, Wds), preserve_range=True, anti_aliasing=True).astype(np.float32)
heM_ds   = resize((he_nuc_can > 0).astype(np.float32),   (Hds, Wds), preserve_range=True).astype(np.float32)
dapiM_ds = resize((dapi_can   > 0).astype(np.float32),   (Hds, Wds), preserve_range=True).astype(np.float32)

# Masked phase correlation (if your skimage supports masks), else unmasked
try:
    shift, err, _ = phase_cross_correlation(dapiF_ds, heF_ds, upsample_factor=100,
                                            reference_mask=dapiM_ds, moving_mask=heM_ds)
except TypeError:
    shift, err, _ = phase_cross_correlation(dapiF_ds, heF_ds, upsample_factor=100)

dy = float(shift[0] * ds)  # +dy moves HE down
dx = float(shift[1] * ds)  # +dx moves HE right
print(f"Translation: dy={dy:.2f}, dx={dx:.2f} (err={err:.6f})")

# ---------------- apply SAME (angle + dx,dy) to nuclei mask and RGB HE, on DAPI canvas ----------------
# 1) H&E nuclei mask (nearest-neighbor)
henuc_rot = rotate(he_nuc, angle=best_angle, resize=True,
                   order=0, mode='constant', cval=0, preserve_range=True).astype(np.float32)
Hr, Wr = henuc_rot.shape
ys = max(0, (Hr - H) // 2); xs = max(0, (Wr - W) // 2)
yd = max(0, (H  - Hr) // 2); xd = max(0, (W  - Wr) // 2)
henuc_can = np.zeros((H, W), dtype=np.float32)
henuc_can[yd:yd+min(H, Hr), xd:xd+min(W, Wr)] = henuc_rot[ys:ys+min(H, Hr), xs:xs+min(W, Wr)]
M = np.array([[1, 0, dx],
              [0, 1, dy]], dtype=np.float32)
henuc_reg = cv2.warpAffine(henuc_can, M, (W, H),
                           flags=cv2.INTER_NEAREST,
                           borderMode=cv2.BORDER_CONSTANT, borderValue=0.0)
henuc_reg = (henuc_reg > 0.5).astype(np.uint8)  # keep binary

# 2) RGB H&E (bilinear)
he_f   = img_as_float32(he)
he_rot = rotate(he_f, angle=best_angle, resize=True,
                order=1, mode='constant', cval=0, preserve_range=True).astype(np.float32)
Hr, Wr = he_rot.shape[:2]
ys = max(0, (Hr - H) // 2); xs = max(0, (Wr - W) // 2)
yd = max(0, (H  - Hr) // 2); xd = max(0, (W  - Wr) // 2)
he_can = np.zeros((H, W, he_rot.shape[2]), dtype=np.float32)
he_can[yd:yd+min(H, Hr), xd:xd+min(W, Wr), :] = he_rot[ys:ys+min(H, Hr), xs:xs+min(W, Wr), :]
he_reg = np.zeros_like(he_can, dtype=np.float32)
for ch in range(he_can.shape[2]):
    he_reg[..., ch] = cv2.warpAffine(he_can[..., ch], M, (W, H),
                                     flags=cv2.INTER_LINEAR,
                                     borderMode=cv2.BORDER_CONSTANT, borderValue=0.0)

# ---------------- save outputs ----------------
imwrite(out_henuc, henuc_reg.astype(np.uint8))
imwrite(out_he, (np.clip(he_reg, 0, 1) * 65535).astype(np.uint16))
print("Wrote:", out_henuc)
print("Wrote:", out_he)

# ---------------- quick QC overlay (DAPI edges in lime over HE) ----------------
plt.figure(figsize=(8,8))
plt.imshow(np.clip(he_reg, 0, 1))
plt.contour(dapi_mask > 0, levels=[0.5], colors=['lime'], linewidths=0.7)
plt.axis('off'); plt.tight_layout()
plt.savefig(out_overlay, dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()
print("Wrote:", out_overlay)