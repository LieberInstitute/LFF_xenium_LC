import os
import glob
import numpy as np
import matplotlib.pyplot as plt

from skimage.io import imread, imsave
from skimage.color import rgb2gray, rgb2lab
from skimage.filters import threshold_otsu
from skimage.filters.rank import entropy
from skimage.morphology import (
    disk,
    remove_small_objects,
    remove_small_holes,
    binary_opening,
    binary_closing
)
from skimage.measure import label, regionprops
from skimage.segmentation import find_boundaries
from skimage.util import img_as_ubyte
from scipy.ndimage import gaussian_filter, binary_fill_holes


# -----------------------------
# Input / output
# -----------------------------
img_dir = "/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing_new/registrations"
out_dir = "/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing_new/NMsegs_adaptive"

os.makedirs(out_dir, exist_ok=True)

files = sorted(glob.glob(os.path.join(img_dir, "Br*_HE.png")))
print("Number of files:", len(files))


# -----------------------------
# Helper functions
# -----------------------------
def normalize01(x, mask=None):
    x = x.astype(np.float32)
    if mask is None:
        lo, hi = np.percentile(x, [1, 99])
    else:
        vals = x[mask]
        lo, hi = np.percentile(vals, [1, 99])
    out = (x - lo) / (hi - lo + 1e-8)
    out = np.clip(out, 0, 1)
    return out


def segment_nm_adaptive(
    f,
    min_obj_size=80,
    min_hole_size=80,
    entropy_disk=5,
    local_bg_sigma=20,
    score_percentile=92
):
    # -----------------------------
    # Read image
    # -----------------------------
img = imread(f)
if img.ndim == 2:
    raw_rgb = np.stack([img, img, img], axis=-1)
else:
    raw_rgb = img[:, :, :3]
if raw_rgb.dtype != np.uint8:
    raw_rgb = raw_rgb.astype(np.float32)
    raw_rgb = raw_rgb / raw_rgb.max()
    raw_rgb = (raw_rgb * 255).astype(np.uint8)
raw_float = raw_rgb.astype(np.float32) / 255.0
R = raw_float[:, :, 0]
G = raw_float[:, :, 1]
B = raw_float[:, :, 2]
gray = rgb2gray(raw_float)
# -----------------------------
# Tissue mask
# -----------------------------
tissue_mask = gray < 0.90
tissue_mask = remove_small_objects(tissue_mask, min_size=100000)
tissue_mask = binary_fill_holes(tissue_mask)
# -----------------------------
# LAB color features
# -----------------------------
lab = rgb2lab(raw_float)
L = lab[:, :, 0]
A = lab[:, :, 1]
BB = lab[:, :, 2]
# Brown/purple-brown chromaticity
lab_brown_score = A + BB
# RGB brownness
rgb_brown_score = R - 0.5 * G - 0.5 * B
# ----------------------------
# Local adaptive darkness
# -----------------------------
local_bg = gaussian_filter(gray, sigma=local_bg_sigma)
relative_dark = local_bg - gray
# -----------------------------
# Texture / granular NM feature
# -----------------------------
gray_u8 = img_as_ubyte(gray)
entropy_img = entropy(gray_u8, disk(entropy_disk))
# -----------------------------
# Normalize features inside tissue
# -----------------------------
lab_brown_norm = normalize01(lab_brown_score, tissue_mask)
rgb_brown_norm = normalize01(rgb_brown_score, tissue_mask)
rel_dark_norm = normalize01(relative_dark, tissue_mask)
entropy_norm = normalize01(entropy_img, tissue_mask)
# -----------------------------
# Combined NM score
# -----------------------------
nm_score = (
    1.2 * lab_brown_norm +
    1.0 * rgb_brown_norm +
    1.0 * rel_dark_norm +
    1.3 * entropy_norm
)
nm_score[~tissue_mask] = 0
# Adaptive threshold
thresh = np.percentile(nm_score[tissue_mask], score_percentile)
nm_mask = nm_score > thresh
# -----------------------------
# Cleanup
# -----------------------------
nm_mask = binary_opening(nm_mask, disk(1))
nm_mask = binary_closing(nm_mask, disk(2))
nm_mask = remove_small_objects(nm_mask, min_size=min_obj_size)
nm_mask = remove_small_holes(nm_mask, area_threshold=min_hole_size)
nm_mask = binary_fill_holes(nm_mask)
# -----------------------------
# Morphology filter
# -----------------------------
lab_mask = label(nm_mask)
clean_mask = np.zeros_like(nm_mask, dtype=bool)
for region in regionprops(lab_mask):
    area = region.area
    solidity = region.solidity
    eccentricity = region.eccentricity
    if area < min_obj_size:
        continue
    # Keep irregular NM bodies and granular regions
    if solidity > 0.30 and eccentricity < 0.98:
        clean_mask[lab_mask == region.label] = True
nm_mask = clean_mask
# -----------------------------
# Overlay
# -----------------------------
boundaries = find_boundaries(nm_mask, mode="outer")
overlay = raw_rgb.copy()
overlay[boundaries] = [0, 255, 0]
return raw_rgb, nm_mask, nm_score, overlay


# -----------------------------
# Run all images
# -----------------------------
for f in files:
    fname = os.path.splitext(os.path.basename(f))[0]
    print("Processing:", fname)

    raw_rgb, nm_mask, nm_score, overlay = segment_nm_adaptive(
        f,
        min_obj_size=80,
        min_hole_size=80,
        entropy_disk=5,
        local_bg_sigma=20,
        score_percentile=92
    )

    imsave(
        os.path.join(out_dir, fname + "_NMmask.png"),
        (nm_mask.astype(np.uint8) * 255)
    )

    imsave(
        os.path.join(out_dir, fname + "_NMoverlay.png"),
        overlay
    )

    plt.figure(figsize=(15, 5))

    plt.subplot(1, 3, 1)
    plt.imshow(raw_rgb)
    plt.title("Raw image")
    plt.axis("off")

    plt.subplot(1, 3, 2)
    plt.imshow(nm_score, cmap="magma")
    plt.title("Adaptive NM score")
    plt.axis("off")

    plt.subplot(1, 3, 3)
    plt.imshow(overlay)
    plt.title("NM boundary overlay")
    plt.axis("off")

    plt.tight_layout()
    plt.savefig(
        os.path.join(out_dir, fname + "_QC.png"),
        dpi=200,
        bbox_inches="tight"
    )
    plt.close()

print("Done.")