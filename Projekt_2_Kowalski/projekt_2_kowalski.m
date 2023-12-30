close all; clear; clc;
a = imread('poj.jpg');
skala = imread('skala_przy_0.jpg');
glaukonit = (a(:,:,1) > 65 &  a(:,:,1) < 80) & (a(:,:,2) > 80 & a(:,:,2) < 100) & (a(:,:,3) > 50 & a(:,:,3) < 68);

SE = strel('disk', 7);
glaukonit = imclose(glaukonit, SE);
glaukonit = medfilt2(glaukonit, [25, 25]);
glaukonit = imfill(glaukonit, "holes");
glaukonit = bwareaopen(glaukonit, 6000);

%subplot(221), imshow(a)
%subplot(222), imshow(imoverlay(a, glaukonit, 'g'));
pole_glaukonit = sum(glaukonit(:)) * (100/186)^2;

b = imread('skrzy_0.jpg');

mika = (b(:,:,1) > 120 & b(:, :, 1) < 230 & b(:,:,2) < 80 & b(:,:,3) > 77 & b(:,:,3) < 220);
SE = strel('disk', 11);
mika = imclose(mika, SE);
mika = imfill(mika, "holes");
mika = bwareaopen(mika, 6000);
pole_mika = sum(mika(:)) * (100/186)^2;

%subplot(121),imshow(b);
%subplot(122),imshow(imoverlay(b, mika, 'r'))

kwarc1 = b(:,:,1) > 8 & b(:,:,1) < 25 & b(:,:,2) > 8 & b(:,:,2) < 25 & b(:,:,3) > 5 & b(:,:,3) < 10;
SE = ones(10);
kwarc1 = imclose(kwarc1, SE);
kwarc1 = imfill(kwarc1, 'holes');
kwarc1 = medfilt2(kwarc1, [25,25]);
kwarc1 = bwareaopen(kwarc1, 12000);

kwarc2 = b(:,:,1) > 170 & b(:,:,2) > 170 & b(:,:,3) > 170;
SE = ones(10);
kwarc2 = imclose(kwarc2, SE);
kwarc2 = imfill(kwarc2, 'holes');
kwarc2 = medfilt2(kwarc2, [25, 25]);
kwarc2 = bwareaopen(kwarc2, 5000);

kwarc3 = b(:,:,1) > 50 & b(:,:,1) < 110 & b(:,:,2) > 50 & b(:,:,2) < 103 & b(:,:,3) > 65 & b(:,:,3) < 95;
SE = ones(10);
kwarc3 = imclose(kwarc3, SE);
kwarc3 = imfill(kwarc3, 'holes');
kwarc3 = medfilt2(kwarc3, [25, 25]);
kwarc3 = bwareaopen(kwarc3, 8000);


kwarc4 = b(:, :, 1) > 120 & b(:,:,1) < 180 & b(:,:,2) > 120 & b(:,:,2) < 150 & b(:,:,3) > 100 & b(:,:,3) < 150;
SE = ones(10);
kwarc4 = imclose(kwarc4, SE);
kwarc4 = imfill(kwarc4, 'holes');
kwarc4 = medfilt2(kwarc4, [25, 25]);
kwarc4 = bwareaopen(kwarc4, 5000);

kwarc5 =  b(:, :, 1) > 50 & b(:,:,1) < 70 & b(:,:,2) > 50 & b(:,:,2) < 70 & b(:,:,3) > 30 & b(:,:,3) < 55;
kwarc5 = imclose(kwarc5, SE);
kwarc5 = imfill(kwarc5, 'holes');
kwarc5 = medfilt2(kwarc5, [25, 25]);
kwarc5 = bwareaopen(kwarc5, 6000);
%subplot(121), imshow(a);
%subplot(122), imshow(imoverlay(a, kwarc, 'y'))

kw = imread('skrzy_60.jpg');
kwa =  kw(:, :, 1) > 150 & kw(:,:,1) < 170 & kw(:,:,2) > 150 & kw(:,:,2) < 170 & kw(:,:,3) > 100 & kw(:,:,3) < 140;
kwa = imclose(kwa, ones(15));
kwa = imfill(kwa, 'holes');
kwa = medfilt2(kwa, [25, 25]);
kwa = bwareaopen(kwa, 10000);
kwa = imrotate(kwa, 60, 'crop');

kwarc = kwarc1 | kwarc2 | kwarc3 | kwarc4 | kwarc5 | kwa;
kwarc = imclose(kwarc, SE);
pole_kwarc= sum(kwarc(:)) * (100/186)^2;

glaukonit = imoverlay(a, glaukonit, 'g');
km = imoverlay(b, kwarc, 'y');
km = imoverlay(km, mika, 'r');

subplot(221), imshow(a), title('Wejściowe zdjęcie 1N')
subplot(222), imshow(glaukonit), title('Znalezione glaukonity - zielone');
subplot(223), imshow(b), title('Wejściowy obraz XN, orientacja 0');
subplot(224), imshow(km), title('Znalezione miki - czerwona i kwarce - żółte');

disp(['Pole glaukonitów: ', num2str(pole_glaukonit), 'µm2'])
disp(['Pole mik: ', num2str(pole_mika), 'µm2'])
disp(['Pole kwarców: ', num2str(pole_kwarc), 'µm2']) 