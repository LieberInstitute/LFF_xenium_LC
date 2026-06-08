import tifffile as tiff
import matplotlib.pyplot as plt

path = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_xenium_LC/processed-data/xenium_imageProcessing/Br0942/xeniumranger_NM_DAPI_Br0942/outs/morphology.ome.tif'

with tiff.TiffFile(path) as tif:
    print(f"Total series: {len(tif.series)}")
    # pick the LAST series (lowest resolution)
    series = tif.series[-1]
    print("axes:", series.axes)
    print("shape:", series.shape)  # should be (11, small_y, small_x)
    img = series.asarray()  # safe now (small!)
# img shape = (Z, Y, X)
Z = img.shape[0]
# plot all planes
fig, axes = plt.subplots(3, 4, figsize=(12, 9))
axes = axes.ravel()
for z in range(Z):
    axes[z].imshow(img[z], cmap='gray')
    axes[z].set_title(f'Z={z}')
    axes[z].axis('off')
# hide extra subplot
for i in range(Z, len(axes)):
    axes[i].axis('off')
plt.tight_layout()
plt.show()