import imopen
import ordfilt2
import imfill
import regionprops
import cv2

image_path = input("Podaj ścieżkę do obrazu: ")
print('1. Regionprops (obraz monochromatyczny)')
print('2. Ordfilt2 (obraz RGB)')
print('3. Ordfilt2 (obraz monochromatyczny)')
print('4. Imopen (obraz logiczny lub monochromatyczny)')
print('5. Imfill (obraz logiczny)')
option = input('Wybierz jedna z opcji (1-5): ')
if option == '1':
    output_path = input('Podaj sciezke do zapisu wynikowego pliku tekstowego: ')
    regionprops.regionprops(image_path, output_path)
    print('Wyniki zostaly zapisane do pliku (z rozszerzeniem).')
elif option == '2':
    rank = int(input('Podaj numer porzadkowy: '))
    mask_size = int(input('Podaj rozmiar maski: '))
    output_path = input('Podaj sciezke do zapisu obrazu: ')
    ordfilt2.ordfilt2_mono(image_path, output_path, rank, mask_size)
    print('Wyniki zostaly zapisane do pliku (z rozszerzeniem).')
elif option == '3':
    rank = int(input('Podaj numer porzadkowy: '))
    mask_size = int(input('Podaj rozmiar maski: '))
    output_path = input('Podaj sciezke do zapisu obrazu: ')
    ordfilt2.ordfilt2_color(image_path, output_path, rank, mask_size)
    print('Wyniki zostaly zapisane do pliku (z rozszerzeniem).')
elif option == '4':
    r = int(input('Podaj rozmiar promienia: '))
    output_path = input('Podaj sciezke do zapisu obrazu: ')
    imopen.imopen(image_path, output_path, r)
    print('Wyniki zostaly zapisane do pliku (z rozszerzeniem).')
elif option == '5':
    output_path = input('Podaj sciezke do zapisu obrazu: ')
    imfill.imfill(image_path, output_path)
    print('Wyniki zostaly zapisane do pliku (z rozszerzeniem).')