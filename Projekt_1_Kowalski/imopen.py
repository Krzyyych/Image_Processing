import numpy as np
import cv2

def create_circle(matrix, r, size):
    center_x = size // 2 - 0.5
    center_y = size // 2 - 0.5

    for i in range(size):
        for j in range(size):
            distance = np.sqrt((i - center_x) ** 2 + (j - center_y) ** 2)
            if distance <= r:
                matrix[i][j] = 1

    return matrix

def erode(image, SE, size):
    center_x = (size - 1) // 2
    center_y = (size - 1) // 2
    erode = image.copy()

    for y in range(image.shape[0]):
        for x in range(image.shape[1]):
            newColor = 255
            for yy in range(-center_y, center_y + 1):
                for xx in range(-center_x, center_x + 1):
                    if SE[yy + center_y][xx + center_x] != 1 or x + xx < 0 or x + xx >= image.shape[1] or y + yy < 0 or y + yy >= image.shape[0]:
                        continue
                    pixel_r = image[y + yy, x + xx]
                    newColor = min(newColor, pixel_r)
            erode[y, x] = newColor

    return erode

def dilate(image, SE, size):
    center_x = (size - 1) // 2
    center_y = (size - 1) // 2
    dilate = image.copy()

    for y in range(image.shape[0]):
        for x in range(image.shape[1]):
            newColor = 0

            for yy in range(-center_y, center_y + 1):
                for xx in range(-center_x, center_x + 1):
                    if SE[yy + center_y][xx + center_x] != 1 or x + xx < 0 or x + xx >= image.shape[1] or y + yy < 0 or y + yy >= image.shape[0]:
                        continue
                    pixel_r = image[y + yy, x + xx]
                    newColor = max(newColor, pixel_r)

            dilate[y, x] = newColor

    return dilate

def imopen(image_path, output_path, r):
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    SE = np.zeros((2 * r, 2 * r), dtype=int)
    SE = create_circle(SE, r, 2 * r)
    result = erode(image, SE, 2 * r)
    result = dilate(result, SE, 2 * r)
    cv2.imwrite(output_path, result)







