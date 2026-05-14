import os
import glob
import numpy as np
from imageio.v3 import imread, imwrite
from PIL import Image

from skimage.color import rgb2gray
from skimage.measure import label, regionprops
from skimage.filters import threshold_otsu, sobel, rank
from skimage.filters.rank import entropy
from skimage.morphology import remove_small_objects, remove_small_holes, binary_opening, binary_closing, disk
from skimage.segmentation import find_boundaries
from scipy.ndimage import binary_fill_holes
import matplotlib.pyplot as plt
from skimage.util import img_as_ubyte

Image.MAX_IMAGE_PIXELS = None

img_dir = "/dcs05/lieber/marmaypag/human_VTA_profiling_LIBD001/human_VTA_spatial/raw-data/images"
out_dir = os.path.join("/dcs05/lieber/marmaypag/human_VTA_profiling_LIBD001/human_VTA_spatial/processed-data/NMsegs", "NM_segmentations_brown")
os.makedirs(out_dir, exist_ok=True)

files = sorted(glob.glob(os.path.join(img_dir, "[0-1][0-9]*.tif")))
print(len(files))

for f in files:
    print(f)
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
# Extract brown NM-like signal
R = raw_float[:, :, 0]
G = raw_float[:, :, 1]
B = raw_float[:, :, 2]
gray = rgb2gray(raw_float)
#gray_8 = (gray * 255).astype(np.uint8)
# Calculate local standard deviation in a small neighborhood
# NM will have high variance; smears will have low variance
#local_std = rank.std(gray, disk(3))

# Normalize and create a mask (threshold may need tuning, e.g., > 10 or 20)
#texture_mask = local_std > 15

# 1) Tissue mask: remove white background
tissue_mask = gray < 0.90
tissue_mask = remove_small_objects(tissue_mask, min_size=1000000)
tissue_mask = binary_fill_holes(tissue_mask)
# 2) NM dark brown/black-ish:
# low intensity + brown enrichment + not hematoxylin-blue/purple
brown_score = (R - B) + 0.5 * (R - G)
brown_score = (brown_score - np.nanmin(brown_score)) / (
    np.nanmax(brown_score) - np.nanmin(brown_score)
)
# Main NM mask: stricter than tissue
brown_mask = (
    tissue_mask &
    (gray < 0.40) &          # dark pixels only; try 0.35-0.55
    (brown_score < 0.6) #&   # brown-ish; try 0.05-0.20
    #(R > B) #&                 # excludes blue/purple nuclei
    #(G > B)
)

#refined_brown_mask = brown_mask & texture_mask
#
#gray_ubyte = img_as_ubyte(gray)
#nm_complexity = entropy(gray_ubyte, disk(5))
#texture_mask = nm_complexity > 5

# Blue / hematoxylin-like nuclei
blue_score = (B - R) + 0.5 * (B - G)
blue_score = (blue_score - np.nanmin(blue_score)) / (
    np.nanmax(blue_score) - np.nanmin(blue_score)
)

blue_mask = (
    tissue_mask &
    (gray < 0.40) &        # nuclei are dark, but can be less dark than NM
    (blue_score > 0.2) &  # blue/purple enrichment
    (B > R)                # stronger blue than red
)

nm_mask = texture_mask & brown_mask | blue_mask

# Cleanup NM only
# remove isolated pixel noise first
nm_mask = binary_opening(nm_mask, disk(2))
# connect nearby NM pixels
nm_mask = binary_closing(nm_mask, disk(1))
# fill within larger NM blobs
nm_mask = binary_fill_holes(nm_mask)
# remove small noisy objects
#nm_mask1 = remove_small_objects(nm_mask, min_size=100)
min_size = 100
max_size = 30000

lab = label(nm_mask)
sizes = np.bincount(lab.ravel())
valid = (sizes >= min_size) & (sizes <= max_size)
valid[0] = False
nm_mask_filtered = valid[lab]

# remove small holes
nm_mask1 = remove_small_holes(nm_mask_filtered, area_threshold=300)



# define region
r1, r2 = 5000, 6000
c1, c2 = 5000, 6000

# Figure 1: mask
plt.figure()
plt.imshow(nm_mask1[r1:r2, c1:c2], cmap='gray')
plt.title("NM mask")
plt.axis('off')

plt.figure()
plt.imshow(texture_mask, cmap='gray')
plt.title("NM mask")
plt.axis('off')
plt.show()
# Figure 2: raw image
plt.figure()
plt.imshow(raw_rgb)
plt.title("Raw image")
plt.axis('off')

plt.show()

boundary = find_boundaries(nm_mask1, mode="outer")
overlay = raw_rgb.copy()
overlay[boundary] = [0, 255, 0]
brain = os.path.splitext(os.path.basename(f))[0]
imwrite(
    os.path.join(out_dir, brain + "_NM_mask.png"),
    (nm_mask1.astype(np.uint8) * 255)
)
overlay_small = overlay[::2, ::2]
imwrite(
    os.path.join(out_dir, brain + "_NM_boundary_overlay.png"),
    overlay_small
)
    