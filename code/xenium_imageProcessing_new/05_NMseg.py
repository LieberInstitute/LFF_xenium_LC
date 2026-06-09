import os
import glob
import numpy as np
import pandas as pd
from imageio.v3 import imread, imwrite
from PIL import Image
import cv2
from skimage.color import rgb2gray
from skimage.measure import label, regionprops_table
from skimage.morphology import (
    remove_small_objects,
    remove_small_holes,
    binary_opening,
    binary_closing,
    disk,
)
from skimage.segmentation import find_boundaries
from scipy.ndimage import binary_fill_holes

Image.MAX_IMAGE_PIXELS = None

img_dir = "/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing_new/registrations"
out_dir = "/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing_new/NMsegs"
os.makedirs(out_dir, exist_ok=True)

files = sorted(glob.glob(os.path.join(img_dir, "Br*_HE.png")))
print("Number of files:", len(files))

for f in files:
    print(f)
    fname = os.path.splitext(os.path.basename(f))[0]
    img = imread(f)
    raw_rgb = img[:,:,:3]
    raw_float = raw_rgb.astype(np.float32) / 255.0
    R = raw_float[:, :, 0]
    G = raw_float[:, :, 1]
    B = raw_float[:, :, 2]
    gray = rgb2gray(raw_float)
    # ----------------------------
    # 1. Tissue mask
    # -----------------------------
    tissue_mask = gray < 0.90
    tissue_mask = remove_small_objects(tissue_mask, min_size=1_000_000)
    tissue_mask = binary_fill_holes(tissue_mask)
    # -----------------------------
    # 2. Brown / NM score
    # -----------------------------
    brown_score = (R - B) + 0.5 * (R - G) 
    brown_score = (brown_score - np.nanmin(brown_score)) / (np.nanmax(brown_score) - np.nanmin(brown_score)+ 1e-8)
    # -----------------------------
    # 3. Strict NM candidate mask
    # -----------------------------
    brown_mask = (
        tissue_mask
        & (gray < 0.25)
        & (brown_score < 0.75)
        & (R > G)
    )
    #brown_mask = (
    # tissue_mask
    # & (gray < 0.35)          # strict darkness filter
    # & (brown_score < 0.8)   # brown enrichment
    # & (R>G)
    # )
    brown_mask = binary_opening(brown_mask, disk(4))
    brown_mask = remove_small_objects(brown_mask, min_size=100)
    brown_mask = binary_fill_holes(brown_mask)
    # -----------------------------
    # 4. Strict blue candidate mask
    # -----------------------------
    # Blue / hematoxylin-like nuclei
    blue_score = (B - R) + 0.5 * (B - G)
    blue_score = (blue_score - np.nanmin(blue_score)) / (
        np.nanmax(blue_score) - np.nanmin(blue_score)
    )
    blue_mask = (
        tissue_mask &
        (gray < 0.30) &        # nuclei are dark, but can be less dark than NM
        (blue_score > 0.5) &  # blue/purple enrichment
        (B > R))
    nm_mask =  brown_mask | blue_mask
    nm_mask = binary_opening(nm_mask, disk(3))
    # -----------------------------
    # 5. Object-level filtering
    # -----------------------------
    lab = label(nm_mask)
    props = regionprops_table(
     lab,
     intensity_image=gray,
     properties=[
         "label",
         "area",
         "eccentricity",
         "solidity",
         "major_axis_length",
         "minor_axis_length",
         "mean_intensity",
     ],
    )
    props = pd.DataFrame(props)
    if len(props) == 0:
        nm_mask1 = np.zeros_like(nm_mask, dtype=bool)
    else:
        keep = (
            (props["area"] >= 100)
            & (props["area"] <= 30000)
            & (props["solidity"] > 0.35))
        keep_labels = props.loc[keep, "label"].values
        nm_mask1 = np.isin(lab, keep_labels)
    # -----------------------------
    # 7. Boundary overlay
    # -----------------------------
    boundary = find_boundaries(nm_mask1, mode="outer")
    overlay = raw_rgb.copy()
    overlay[boundary] = [0, 255, 0]
    brain = os.path.basename(f).split('_HE')[0]
    imwrite(os.path.join(out_dir, brain + "_NM_mask.png"),(nm_mask1.astype(np.uint8) * 255))
    new_w = int(img.shape[1] * 0.3)
    new_h = int(img.shape[0] * 0.3)
    overlay_small = cv2.resize(overlay, (new_w, new_h), interpolation=cv2.INTER_AREA)
    imwrite(os.path.join(out_dir, brain + "_NM_boundary_overlay.png"),overlay)


#import matplotlib.pyplot as plt
#
#
#boundary = find_boundaries(nm_mask1, mode="outer")
#overlay = raw_rgb.copy()
#overlay[boundary] = [0, 255, 0]
#
#r1, r2 = 10400, 11500
#c1, c2 = 3500, 4700
## Figure 1: mask
#plt.figure()
#plt.imshow(overlay[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)
#
#r3, r4 = 23500, 23800
#c3, c4 = 1800, 2100
## Figure 1: mask
#plt.figure()
#plt.imshow(overlay[r3:r4, c3:c4], cmap='gray')
#plt.show(block=False)
#
#
#plt.figure()
#plt.imshow(overlay, cmap='gray')
#plt.show(block=False)
#
#plt.figure()
#plt.imshow(brown_mask[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)
#
#plt.figure()
#plt.imshow(gray_adapt[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)
#
#plt.figure()
#plt.imshow(nm_mask1[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)
#
#plt.figure()
#plt.imshow(img[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)

#r1, r2 = 12000, 13000
#c1, c2 = 8000, 9200
#
#plt.figure()
#plt.imshow(brown_mask[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)
#
#plt.figure()
#plt.imshow(blue_mask[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)
#
#plt.figure()
#plt.imshow(img[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)
#
#plt.figure()
#plt.imshow(nm_mask1[r1:r2, c1:c2], cmap='gray')
#plt.show(block=False)