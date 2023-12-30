import numpy as np
import cv2

def dilation(image, mask):
    result = np.zeros_like(image)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            result[i, j] = max(image, i, j, mask, 1)
    return result

def max(img, i, j, mask, radius):
    for i0 in range(i - radius, i + radius + 1):
        for j0 in range(j - radius, j + radius + 1):
            if not (i0 < 0 or j0 < 0 or i0 >= img.shape[0] or j0 >= img.shape[1]):
                if ((i0 - i) * 2 + (j0 - j) * 2) <= radius ** 2:
                    if mask[i0, j0] == 1:
                        return 1
    return 0

def imfill(image_path, output_path):
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    maska = np.zeros(image.shape)
    for i in range(image.shape[0]):
        for j in range(image.shape[1]):
            if i == 0 or j == 0 or i == image.shape[0] - 1 or j == image.shape[1] - 1:
                maska[i, j] = 1

    neg = np.logical_not(image)

    while True:
        temp = dilation(image, maska) & neg

        if np.array_equiv(temp, maska):
            cv2.imwrite(output_path, np.logical_not(maska) * 255)
            return
        x = 0
        maska = temp


