from aicsimageio import AICSImage
from skimage.transform import resize
from skimage import io 
import numpy as np

# Path to the input and output files
input_path = "/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/exp66/2024_03_20_Evelien_BMK_mouse_brain_Axioscan7/BMK_experiment_20240320-1/BMK_experiment_20240320-1.tif"
output_path = "/data/projects/s00/2024-07-UTBenchmark-SpC/data/raw/exp66/BMK_experiment_20240320-1_downsampled.tif"

# Read the image using AICSImageIO
img = AICSImage(input_path)
original_image = img.data

# Assuming a 2D image; if multi-channel, you might need to select a specific channel or merge channels
image_data = original_image[0, 0, 0, :, :, :]  # Adjust indexing as needed based on image dimensions

# Determine downsample factor
downsample_factor = 0.65
new_size = (int(image_data.shape[0] * downsample_factor), int(image_data.shape[1] * downsample_factor), 3)

# Downsample the image
downsampled_image = resize(image_data, new_size, anti_aliasing=True, preserve_range=True).astype(np.uint8)

# Save the downsampled image using skimage
io.imsave(output_path, downsampled_image)

print(f"Downsampled image saved to: {output_path}")
