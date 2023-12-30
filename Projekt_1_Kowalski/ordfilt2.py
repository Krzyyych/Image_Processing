import numpy as np
import cv2

def ordfilt2_mono(image_path, output_path, rank, mask_size):
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    height, width = image.shape[:2]
    margin = mask_size // 2
    padded_image = np.pad(image, margin, mode='constant')
    result = np.zeros_like(image)

    for y in range(height):
        for x in range(width):
            patch = padded_image[y:y+mask_size, x:x+mask_size]
            sorted_patch = np.sort(patch.flatten())
            result[y, x] = sorted_patch[rank-1]
    cv2.imwrite(output_path, result)

def ordfilt2_color(image_path, output_path, rank, mask_size):
    image = cv2.imread(image_path)
    height, width = image.shape[:2]
    num_channels = image.shape[2]
    margin = mask_size // 2
    padded_image = np.pad(image, ((margin, margin), (margin, margin), (0, 0)), mode='constant')
    result = np.zeros_like(image)

    for c in range(num_channels):
        for y in range(height):
            for x in range(width):
                patch = padded_image[y:y+mask_size, x:x+mask_size, c]
                sorted_patch = np.sort(patch.reshape(mask_size*mask_size))
                result[y, x, c] = sorted_patch[rank-1]

    cv2.imwrite(output_path, result)




