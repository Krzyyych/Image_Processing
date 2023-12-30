import cv2
import math

def regionprops(image_path, output_path):
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)

    max = 255
    centroids = [{'X': 0, 'Y': 0, 'Count': 0} for i in range(max)]
    equivalent_diameters = [0] * max

    for i in range(max):
        centroids[i]['X'] = 0
        centroids[i]['Y'] = 0
        centroids[i]['Count'] = 0

    for x in range(image.shape[1]):
        for y in range(image.shape[0]):
            id = image[y, x]
            centroids[id]['X'] += x + 1
            centroids[id]['Y'] += y + 1
            centroids[id]['Count'] += 1

    for i in range(max):
        if centroids[i]['Count'] == 0:
            centroids[i]['X'] = -1
            centroids[i]['Y'] = -1
            equivalent_diameters[i] = 0
        else:
            centroids[i]['X'] /= centroids[i]['Count']
            centroids[i]['Y'] /= centroids[i]['Count']
            equivalent_diameters[i] = math.sqrt(4 * centroids[i]['Count'] / math.pi)

    with open(output_path, 'w') as f:
        f.write('Nr,Centroids,EquivalentDiameters\n')
        for i in range(max):
            f.write(f'{i},{centroids[i]["X"]:.2f},{centroids[i]["Y"]:.2f},{equivalent_diameters[i]:.2f}\n')

