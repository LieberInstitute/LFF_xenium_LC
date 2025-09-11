import tifffile as tf
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

# ---- INPUTS ----
img_path = '/dcs05/lieber/marmaypag/LFF_spatialLC_LIBD4140/LFF_spatial_LC/processed-data/xenium_imageProcessing/Br6538_HE.png'
img = np.array(Image.open(img_path))
ome_tif_path = img_path.replace('.png', '.ome.tif')

pixel_size_um = 0.25  # microns per pixel

# ---- METADATA ----
metadata = {
    'PhysicalSizeX': pixel_size_um,
    'PhysicalSizeXUnit': 'µm',
    'PhysicalSizeY': pixel_size_um,
    'PhysicalSizeYUnit': 'µm'
}

# Compression options
options = dict(
    photometric='rgb',
    tile=(1024, 1024),
    compression='jpeg',
    resolution=(1e4 / pixel_size_um, 1e4 / pixel_size_um),
    resolutionunit='CENTIMETER',
    metadata=metadata
)

# Write base + pyramid levels
with tf.TiffWriter(ome_tif_path, bigtiff=True) as tif:
    # Write level 0
    tif.write(img, subifds=4, **options)
    # Write 4 pyramid levels
    for i in range(4):
        img = img[::2, ::2]  # downsample by 2
        tif.write(img, subfiletype=1, **options)

print(f"Saved: {ome_tif_path}")