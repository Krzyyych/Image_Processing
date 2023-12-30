%% Laboratoria 1 - podstawowe działanie na obrazach szarych
close all; clear; clc;

%obrazy w MatLabie są zapisywany jako inty (inaczej niż wszystko), ze
%względu na oszczędność miejsca
%obrazem moze byc nie tylko zdjęcie ale też np. rozkład gęstości
%temperatury

a = imread('cameraman.tif');
b = a - 50;     %odjete od wszystkich pikseli powoduje przyciemnienie
c = a + 50;     %dodanie do wszystkich pikseli powoduje rozjaśnienie
d = a * 3;      %pomnożenie wszystkich pikseli powoduje skrajne rozjąsnienie
subplot(221), imshow(a);
subplot(222), imshow(b);
subplot(223), imshow(c);
subplot(224), imshow(d);

figure;
e = double(a)/255;      %przejscie z obrazu w intach na double
f = uint8(e*255);       %przejscie z obrazu w doublach na inty
g = uint8(e) * 255;     %powoduje ze otrzymujemy obraz czarno bialy
subplot(221), imshow(e);
subplot(222), imshow(f);
subplot(223), imshow(g);

%obraz w intach jest z zakresu 0 - 255
%obraz w doublach jest z zakresu 0 - 1
%wyjscie ponizej zakres w matlabie powoduje zamiane na 0
%wyjscie ponizej zakresu w C++ powoduje przeskoczenie z czarnego na bialy

%% Laboratoria 1 - podstawowe działania na obrazach kolorowych
close all; clear; clc;

%obraz kolorowy rgb, sklada sie z 3 wymiarowego uint8, kolory czerwony,
%zielony i niebieski
a = imread('onion.png');
subplot(221), imshow(a);
for k = 1 : 3
    subplot(2, 2, k + 1), imshow(a(:,:,k)); %wyswietlenie wszystkich skladowych rgb oddzielnie
end

%obraz kolorowy moze byc zapisany jako obraz indeksowany (opłacalne gdy
%mamy mało kolorów)
[map, leg] = rgb2ind(a, 1200);  %drugi argument to ilosc kolorow
b = ind2rgb(map, leg);  %zamiast leg mozna od razu podac jet (palete)
figure;
subplot(121), imshow(a);
subplot(122), imshow(b);

%% Laboratoria 1 - analiza kolorow z imtool
close all; clear; clc;
a = imread('cameraman.tif');
%w konsoli odpalamy imtool
%do analizy obrazow kolorwych: 
%druga ikona z lewej do znajdowania wartosci koloru
%linijka do mierzenia odleglosci (wiemy ze np. 1mm oznacza 43 piksele)

%do analizy obrazow czarnobialych:
%ostatnia opcja choose colormap i wybieramy interesujaca nas
%moze to pomoc uwydatnic elementy lub dla konwencji (niska temperatura
%kolor niebieski, a wysoko kolor czerwony)

b = imfinfo('cameraman.tif'); %wyswietla informacje o obrazie
c = regionprops(a, 'all'); %etykietowanie 
%drugi argument "all" mozna zmienic na "area" - tylko interesujacy nas
%centroit - liczy tylko dla srodka masy
%boundingbox - najmniejszy prostokat opisany na danym ksztalcie
%dla ulatwienia zmiana figury wkleslej na wypukla

%% Laboratoria 1 - przekroj przez obraz kolorow rgb
close all; clear; clc;
a = imread('onion.png');
[Nz, Nx, k] = size(a); %musimy uwzglednic k dla obrazow kolorowych
subplot(121), imshow(a);
subplot(122), improfile(a, [1, Nx], [1, Nz]); %wykres wartosci, przekroju wartosci rgb

%% Laboratoria 1 - przekroj przez obraz kolorow rgb II
close all; clear; clc;
a = imread('onion.png');
[Nz, Nx, k] = size(a); %musimy uwzglednic k dla obrazow kolorowych
subplot(121), imshow(a);
b = improfile(a, [1, Nx/2, Nx, Nx, 1], [1, Nz/2, 1, Nz, Nz]);
%subplot(122), plot(b(:, :, 2), 'g');
%hold on 
%plot(b(:, :, 3), 'b');
N = size(b, 1);
subplot(122), plot(1:N, b(:, 1, 2), 'g', 1:N, b(:, 1, 3), 'b');

%% Laboratoria 1 - metody interpolacji do zmiany rozmiaru obrazu
%rodzaje interpolacji:
%najblizszego sasiada
%bilinowa
%bikubiczna

close all; clear; clc;
a = checkerboard(8, 4, 4); %szachownica
a = imread('cameraman.tif');
skala = 0.7;
b = imresize(a, skala, 'nearest');
c = imresize(a, skala, 'bilinear');
d = imresize(a, skala,'bicubic');
subplot(221), imshow(a);
subplot(222), imshow(b);
subplot(223), imshow(c);
subplot(224), imshow(d);

%% Laboratoria 2 - korekta gamma
close all; clear; clc;
a = imread("cameraman.tif");
gamma = [1/8, 1/4, 1/2, 1, 2, 4];
a = double(a)/255;
for k = 1:6
    b = a.^gamma(k);
    subplot(2,3,k);
    imshow(b);
    title(["\gamma = ", num2str(gamma(k))]);
end

%gamma mniejsza od 1 wyroznia nam szczegoly z ciemnych fragmentow obrazu
%gamma wieksza od 1 wyroznia nam szczegoly z jasnych fragmentow obrazu
%korekta gamma jest jednym z podstawowych wspolczynnikow w monitorach

%% Laboratoria 2 - normalizacja
%czyli bedziemy rozscigac kolory ale zostawiac proporcje
close all; clear; clc;
a = imread('pout.tif');
subplot(221), imshow(a);
subplot(222), imhist(a, 255);
b = imadjust(a);
subplot(223), imshow(b);
subplot(224), imhist(b, 255);

%% Laboratoria 2 - wyrownanie histogramu
%liczymy kumulante histogramu i normalizujemy ja do przedzialu 0 - 1
%dzielimy os pionowa (y) na jakas liczbe wartosci np. 4
%rzutujemy te 4 wartosci na kumulante oraz os pozioma (x)
%wszystkie piksele z kolejnych przedzialow maja 0, 1/3, 2/3, 1
close all; clear; clc;
a = imread('pout.tif');
subplot(221), imshow(a);
subplot(222), imhist(a, 256);
b = histeq(a, 16);
subplot(223), imshow(b);
subplot(224), imhist(b, 256);

%% Laboratoria 2 - Clahe adaptacyjne wyrownanie histogramu
%podobnie jak wyzej, ale komulanta tylko dla naszego wycinka
%min i max tak samo tylko ze z wycinka
close all; clear; clc;
a = imread('saturn.png');
a = rgb2gray(a);
subplot(221), imshow(a);
subplot(222), imhist(a, 256);
b = adapthisteq(a, 'ClipLimit', 0.05);
subplot(223), imshow(b);
subplot(224), imhist(b, 256);

%% Laboratoria 2 - binearyzacja przypadek prostszy
%zamiana obrazu na obraz logiczny true, false
close all; clear; clc;
a = imread('coins.png');
subplot(131), imshow(a);
b = a > 90;
subplot(132), imshow(b);
bin = medfilt2(b, [3 3]);
subplot(133), imshow(bin);

%% Laboratoria 2 - binearyzacja przypadek trudniejszy
close all; clear; clc;
a = imread('peppers.png');
subplot(221), imshow(a);
r = a(:,:,1);
g = a(:,:,2);
b = a(:,:,3);
binr = r > 160;
binb = b > 140;
bing = g > 110;
bin = binr & binb & bing;
subplot(222), imshow(bin);
bin = bwareaopen(bin, 500); %usuniecie malych obiektow, o mniejszej powierzchni niz drugi argument
subplot(223), imshow(bin);
bin = medfilt2(bin, [3 3]);
subplot(224), imshow(bin);

%% Laboratoria 2 - geometryczne
close all; clear; clc;
a = imread('cameraman.tif');
b = circshift(a, [50, 100]);               %przesuniecie
c = imrotate(a, 30, 'loose', 'bilinear');  %rotacja (obraz, kat(w stopniach), 'loose' - rozmiar obrazka sie zwiekszy, interpolacja)
d = imrotate(a, 30, 'crop', 'bilinear');   %rotacja (obraz, kat(w stopniach), 'crop' - jesli chcemy przyciac do tego samego rozmiaru, interpolacja)
e = flipud(a);                             %obrocenie gora dol; 
f = fliplr(a);                             %obrocenie lewo prawo
g = padarray(a , [50, 100], 'replicate', 'post'); %doklejenie kolumn na poczatku i koncu, wejsciowo 12345, replicate - 111111234555555, 
h = padarray(a , [50, 100], 'symmetric', 'post'); %z 12345 robi 543211234554321
i = padarray(a, [50, 100], 'circular', 'post'); %z 12345 circular 123451234512345
subplot(331), imshow(a);
subplot(332), imshow(b);
subplot(333), imshow(c);
subplot(334), imshow(d);
subplot(335), imshow(e);
subplot(336), imshow(f);
subplot(337), imshow(g);
subplot(338), imshow(h);
subplot(339), imshow(i);

%% Laboratoria 2 - przeksztalcenie afiniczne
%macierz przeksztalcen musi miec wyznacznik rozny od 0
% [a b 0]
% [c d 0]
% [0 0 1]
% a - w poziomie
% b - w pionie
% c - ukosne przesuniecie kopniecie
% d - ukosne przesuniecie kopniecie
close all; clear; clc;
a = imread('cameraman.tif');
mac1 = affine2d([1 2 0; 2 1 0; 0 0 1]);
mac2  = projective2d([1 0 0; 0 1 0; 0 0 6]);
b = imwarp(a, mac1);
c = imwarp(a, mac2);
subplot(131), imshow(a); 
subplot(132), imshow(b);
subplot(133), imshow(c);

%% Laboratoria 3 - filtracje liniowe dolnoprzepustowe jednorodna uśredniajaca
%ostre krawedzi ulegaja wyplaszczeniu
%przez efekt brzegowy pojawia sie czarna ramka (przez uzupelnianie)
%domyslnie jest uzupelniane zerami, ale mozna zamienic np symmetric
close all; clear; clc;
a = imread('cameraman.tif');
N = 7;
maska1 = ones(N) / (N * N);
maska2 = ones(N + 10) / ((N + 10) * (N + 10));
b = imfilter(a, maska1);
c = imfilter(a, maska1, 'symmetric');
d = imfilter(a, maska2, "symmetric");
subplot(221), imshow(a);
subplot(222), imshow(b);
subplot(223), imshow(c);
subplot(224), imshow(d);

%% Laboratoria 3 - filtracja liniowa dolnoprzepustowa uśredniajaca ważona maską gaussa
%zbyt duze maski sprowadzaja sie do delty diraca
close all; clear; clc;
a = imread('cameraman.tif');
N = 7;
maska = fspecial('gaussian', [1 1], N/4);
b = imfilter(a, maska, 'symmetric');
subplot(121), imshow(a);
subplot(122), imshow(b);

%% Laboratoria 3 - filtracja liniowa górnoprzepustowa maska prewitta
close all; clear; clc
a = imread('cameraman.tif');
a = double(a)/255;
maska = [1 1 1; 0 0 0; -1 -1 -1];
b = abs(imfilter(a, maska, 'symmetric'));
c = abs(imfilter(a, maska', 'symmetric'));
d = sqrt(b.^2 + c.^2); %magnituda krawędzi
subplot(221), imshow(a);
subplot(222), imshow(b);
subplot(223), imshow(c);
subplot(224), imshow(d);

%% Laboratoria 3 - filtracja liniowa górnoprzepustowa maska sobela
close all; clear; clc;
a = imread('cameraman.tif');
a = double(a)/255;
maska = [1 2 1; 0 0 0; -1 -2 -1];
b = abs(imfilter(a, maska, 'symmetric'));
c = abs(imfilter(a, maska', 'symmetric'));
d = sqrt(b.^2 + c.^2); %magnituda krawędzi
subplot(221), imshow(a);
subplot(222), imshow(b);
subplot(223), imshow(c);
subplot(224), imshow(d);

%% Laboratoria 3 - zadanko z maskami
close all; clear; clc
a = zeros(120);
a(21:100, 21:100) = 1;
maska = [1 0 -1];
b = abs(imfilter(a, maska) .* imfilter(a, maska'));
b = b & a;
c = uint8(a + b);
leg = [0 0 0; 1 1 1; 0 0 1];
imshow(c, leg);

%% Laboratoria 3 - filtracja liniowa górnoprzepustowa maska unsharp
% parametr alfa do sterowania moca
close all; clear; clc;
a = imread('peppers.png');
a = double(a)/255;
maska1 = [0 -1 0; -1 4 -1; 0 -1 0];
maska2 = [0 -1 0; -1 5 -1; 0 -1 0];
b = imfilter(a, maska1, 'symmetric');
c = imfilter(a, maska2, "symmetric");
d = sqrt(b(:,:,1).^2 + b(:,:,2).^2 + b(:,:,3).^2);
subplot(221), imshow(a);
subplot(222), imshow(b);
subplot(223), imshow(c);
subplot(224), imshow(d);

%% Laboratoria 3 - filtracja nieliniowa medianowa
close all; clear; clc;
a = imread('cameraman.tif');
b = medfilt2(a, [5 5], 'symmetric');
subplot(121), imshow(a);
subplot(122), imshow(b);

%% Laboratoria 3 - filtracja nieliniowa adaptacyjna wienera
close all; clear; clc;
a = imread('cameraman.tif');
b = wiener2(a, [5 5]);
subplot(121), imshow(a);
subplot(122), imshow(b);

%% Laboratoria 3 - zaszumianie i odszumianie
close all; clear; clc;
a = imread('cameraman.tif');
as = imnoise(a, 'gaussian');
N = 5;
jednorodny = imfilter(as, ones(N) / (N * N), 'symmetric');
gauss_maska = fspecial('gaussian', [1 1], N/4);
gauss = imfilter(as, gauss_maska, 'symmetric');
median = medfilt2(as, [5 5], 'symmetric');
wiener = wiener2(as, [5 5]);
subplot(231), imshow(a);
title('Zwykly cameraman');

subplot(232), imshow(as);
title("Cameraman zaszumiony");

subplot(233), imshow(jednorodny);
title("Odszumiony jednorodnym");

subplot(234), imshow(gauss);
title("Odszumiony gaussem");

subplot(235), imshow(median);
title("Odszumiony medianowym");

subplot(236), imshow(wiener);
title("Odszumiony wienerem");

%% Laboratoria 3 - szum impulsowy
close all; clear; clc;
a = imread('cameraman.tif');
as = imnoise(a, 'salt & pepper');
N = 5;
jednorodny = imfilter(as, ones(N) / (N * N), 'symmetric');
gauss_maska = fspecial('gaussian', [1 1], N/4);
gauss = imfilter(as, gauss_maska, 'symmetric');
median = medfilt2(as, [3 3], 'symmetric');
wiener = wiener2(as, [5 5]);
subplot(231), imshow(a);
title('Zwykly cameraman');

subplot(232), imshow(as);
title("Cameraman zaszumiony");

subplot(233), imshow(jednorodny);
title("Odszumiony jednorodnym");

subplot(234), imshow(gauss);
title("Odszumiony gaussem");

subplot(235), imshow(median);
title("Odszumiony medianowym");

subplot(236), imshow(wiener);
title("Odszumiony wienerem");

%% Laboratoria 3 - szum poissona tylko dla int
close all; clear; clc;
a = imread('cameraman.tif');
as = imnoise(a, 'poisson');
N = 5;
jednorodny = imfilter(as, ones(N) / (N * N), "symmetric");
gauss_maska = fspecial('gaussian', [1 1], N/4);
gauss = imfilter(as, gauss_maska, 'symmetric');
median = medfilt2(as, [3 3], 'symmetric');
wiener = wiener2(as, [5 5]);
subplot(231), imshow(a);
title('Zwykly cameraman');

subplot(232), imshow(as);
title("Cameraman zaszumiony");

subplot(233), imshow(jednorodny);
title("Odszumiony jednorodnym");

subplot(234), imshow(gauss);
title("Odszumiony gaussem");

subplot(235), imshow(median);
title("Odszumiony medianowym");

subplot(236), imshow(wiener);
title("Odszumiony wienerem");

%% Laboratoria 3 - szum speckle
close all; clear; clc;
a = imread('cameraman.tif');
as = imnoise(a, "speckle");
N = 5;
jednorodny = imfilter(as, ones(N) / (N * N), 'symmetric');
gauss_maska = fspecial('gaussian', [1 1], N/4);
gauss = imfilter(as, gauss_maska, 'symmetric');
median = medfilt2(as, [3 3], 'symmetric');
wiener = wiener2(as, [5 5]);
subplot(231), imshow(a);
title('Zwykly cameraman');

subplot(232), imshow(as);
title("Cameraman zaszumiony");

subplot(233), imshow(jednorodny);
title("Odszumiony jednorodnym");

subplot(234), imshow(gauss);
title("Odszumiony gaussem");

subplot(235), imshow(median);
title("Odszumiony medianowym");

subplot(236), imshow(wiener);
title("Odszumiony wienerem");

%% Laboratoria 3 - VMF i AMF
% VMF tylko dla obrazow kolorowych, bedziemy robic w c++
close all; clear; clc;
a = imread('cameraman.tif');
b = edge(a, 'canny');   %wyszukiwanie krawedzi canny
c = edge(a, 'prewitt'); %wyszukiwanie krawedzi prewitt
d = ordfilt2(a, 1, ones(7)); %filtr porządkujący
e = stdfilt(a, ones(7)); %filtr odchylenia standardowego
f = entropyfilt(a, ones(9));
subplot(231), imshow(a);
title("Zwykly cameraman");

subplot(232), imshow(b);
title("Canny");

subplot(233), imshow(c);
title("Prewitt");

subplot(234), imshow(d);
title("Porzadkujacy");

subplot(235), imshow(e);
title("Odchylenie standardowe");

subplot(236), imagesc(f); axis image; colorbar('vertical');
title('Entropy');

%% Laboratoria 3 - symulacja ruchu podczas robienia zdjecia
close all; clear; clc;
a = imread('cameraman.tif');
maska = fspecial('motion', 11, -30);
b = imfilter(a, maska, 'symmetric');
c = deconvblind(b, maska);
d = deconvreg(b, maska);
e = deconvwnr(b, maska);
f = deconvlucy(b, maska);
subplot(231), imshow(a);
title("Zwykly cameraman");

subplot(232), imshow(b);
title("Ruch");

subplot(233), imshow(c);
title("Deconvblind");

subplot(234), imshow(d);
title("Deconvreg");

subplot(235), imshow(e);
title("Deconvwnr");

subplot(236), imshow(f);
title('Deconvlucy');

%% Laboratoria 4 - erozja
%erozja jest addytywna
close all; clear; clc;
SE1 = strel('disk', 5);
SE2 = strel('line', 10, 30);
SE3 = strel('arbitrary', [0 1 0; 1 1 1; 0 1 0]);
a = imread('cameraman.tif');
b = imread('circles.png');
c = imread('peppers.png');
e1 = imerode(a, SE1);
e2 = imerode(a, SE2);
e3 = imerode(a, SE3);
e4 = imerode(b, SE1);
e5 = imerode(b, SE2);
e6 = imerode(b, SE3);
e7 = imerode(c, SE1);
e8 = imerode(c, SE2);
e9 = imerode(c, SE3);

subplot(221), imshow(a);
title("Cameraman");
subplot(222), imshow(e1);
title('Element jako dysk');
subplot(223), imshow(e2);
title('Element jako linia');
subplot(224), imshow(e3);
title('Element jako arbitary');

figure;
subplot(221), imshow(b);
title("Circles");
subplot(222), imshow(e4);
title('Element jako dysk');
subplot(223), imshow(e5);
title('Element jako linia');
subplot(224), imshow(e6);
title('Element jako arbitary');

figure;
subplot(221), imshow(c);
title("Papper");
subplot(222), imshow(e7);
title('Element jako dysk');
subplot(223), imshow(e8);
title('Element jako linia');
subplot(224), imshow(e9);
title('Element jako arbitary');

%% Laboratoria 4 - zadanko z erozją
% dla circle.png erozja elementem SE1 = ones(11)
% dla tego samego obrazka seria erozji elementem SE2 = ones(3)
% ile iteracji elementem SE2 jest potrzebnych aby osiagnac to samo co SE1?
close all; clear; clc;
a = imread('circles.png');
SE1 = ones(11);
SE2 = ones(3);
e1 = imerode(a, SE1);
N = 1;
e2 = imerode(a, SE2);
while ~isequal(e1, e2)
    N = N + 1;
    e2 = imerode(e2, SE2);
end
disp(N);

%% Laboratoria 4 - dylacja/dylatacja
%dylacja jest addytywna
close all; clear; clc;
SE1 = strel('disk', 5);
SE2 = strel('line', 10, 30);
SE3 = strel('arbitrary', [0 1 0; 1 1 1; 0 1 0]);
a = imread('cameraman.tif');
b = imread('circles.png');
c = imread('peppers.png');
e1 = imdilate(a, SE1);
e2 = imdilate(a, SE2);
e3 = imdilate(a, SE3);
e4 = imdilate(b, SE1);
e5 = imdilate(b, SE2);
e6 = imdilate(b, SE3);
e7 = imdilate(c, SE1);
e8 = imdilate(c, SE2);
e9 = imdilate(c, SE3);

subplot(221), imshow(a);
title("Cameraman");
subplot(222), imshow(e1);
title('Element jako dysk');
subplot(223), imshow(e2);
title('Element jako linia');
subplot(224), imshow(e3);
title('Element jako arbitary');

figure;
subplot(221), imshow(b);
title("Circles");
subplot(222), imshow(e4);
title('Element jako dysk');
subplot(223), imshow(e5);
title('Element jako linia');
subplot(224), imshow(e6);
title('Element jako arbitary');

figure;
subplot(221), imshow(c);
title("Papper");
subplot(222), imshow(e7);
title('Element jako dysk');
subplot(223), imshow(e8);
title('Element jako linia');
subplot(224), imshow(e9);
title('Element jako arbitary');

%% Laboratoria 4 - otwarcie morfologiczne
% otwarcie jest dylacja erozji tym samym elementem
close all; clear; clc;
SE1 = strel('disk', 5);
SE2 = strel('line', 10, 30);
SE3 = strel('arbitrary', [0 1 0; 1 1 1; 0 1 0]);
a = imread('cameraman.tif');
b = imread('circles.png');
c = imread('peppers.png');
e1 = imopen(a, SE1);
e2 = imopen(a, SE2);
e3 = imopen(a, SE3);
e4 = imopen(b, SE1);
e5 = imopen(b, SE2);
e6 = imopen(b, SE3);
e7 = imopen(c, SE1);
e8 = imopen(c, SE2);
e9 = imopen(c, SE3);

subplot(221), imshow(a);
title("Cameraman");
subplot(222), imshow(e1);
title('Element jako dysk');
subplot(223), imshow(e2);
title('Element jako linia');
subplot(224), imshow(e3);
title('Element jako arbitary');

figure;
subplot(221), imshow(b);
title("Circles");
subplot(222), imshow(e4);
title('Element jako dysk');
subplot(223), imshow(e5);
title('Element jako linia');
subplot(224), imshow(e6);
title('Element jako arbitary');

figure;
subplot(221), imshow(c);
title("Papper");
subplot(222), imshow(e7);
title('Element jako dysk');
subplot(223), imshow(e8);
title('Element jako linia');
subplot(224), imshow(e9);
title('Element jako arbitary');

%% Laboratoria 4 - zadanko z otwarciem
% z obrazu blobs.png usunac krawędzie pionowe
% poziome mają zostać, a ukośne jak wyjdzie
close all; clear; clc;
a = imread('blobs.png');
SE1 = strel('line', 10, 0);
SE2 = ones(1, 8);
o1 = imopen(a, SE1);
o2 = imopen(a, SE2);
subplot(131), imshow(a);
subplot(132), imshow(o1);
subplot(133), imshow(o2);

%% Laboratoria 4 - zamknięcie morfologiczne
% zamkniecie jest erozja dylacji tym samym elementem
% erozja <= otwarcie  | zamkniecie <= dylacja
close all; clear; clc;
SE1 = strel('disk', 5);
SE2 = strel('line', 10, 30);
SE3 = strel('arbitrary', [0 1 0; 1 1 1; 0 1 0]);
a = imread('cameraman.tif');
b = imread('circles.png');
c = imread('peppers.png');
e1 = imclose(a, SE1);
e2 = imclose(a, SE2);
e3 = imclose(a, SE3);
e4 = imclose(b, SE1);
e5 = imclose(b, SE2);
e6 = imclose(b, SE3);
e7 = imclose(c, SE1);
e8 = imclose(c, SE2);
e9 = imclose(c, SE3);

subplot(221), imshow(a);
title("Cameraman");
subplot(222), imshow(e1);
title('Element jako dysk');
subplot(223), imshow(e2);
title('Element jako linia');
subplot(224), imshow(e3);
title('Element jako arbitary');

figure;
subplot(221), imshow(b);
title("Circles");
subplot(222), imshow(e4);
title('Element jako dysk');
subplot(223), imshow(e5);
title('Element jako linia');
subplot(224), imshow(e6);
title('Element jako arbitary');

figure;
subplot(221), imshow(c);
title("Papper");
subplot(222), imshow(e7);
title('Element jako dysk');
subplot(223), imshow(e8);
title('Element jako linia');
subplot(224), imshow(e9);
title('Element jako arbitary');

%% Laboratoria 4 - gradient morfologiczny
% są 4 rodzaje gradientu
% 1. dylacja - obraz
% 2. obraz - erozja
% 3. dylacja - erozja
% 4. (dylacja - erozja) / 2
close all; clear; clc;
a = imread('cameraman.tif');
b = imread('peppers.png');
c = imread('circles.png');
SE = strel('disk', 5);
a1 = imdilate(a, SE) - a;
a2 = a - imerode(a, SE);
a3 = imdilate(a, SE) - imerode(a, SE);
a4 = (imdilate(a, SE) - imerode(a, SE))/2;
b1 = imdilate(b, SE) - b;
b2 = b - imerode(b, SE);
b3 = imdilate(b, SE) - imerode(b, SE);
b4 = (imdilate(b, SE) - imerode(b, SE))/2;
c1 = imdilate(c, SE) - c;
c2 = c - imerode(c, SE);
c3 = imdilate(c, SE) - imerode(c, SE);
c4 = (imdilate(c, SE) - imerode(c, SE))/2;
subplot(231), imshow(a);
title('Cameraman');
subplot(232), imshow(a1);
title('Dylacja - Obraz');
subplot(233), imshow(a2);
title('Obraz - Erozja');
subplot(234), imshow(a3);
title('Dylacja - Erozja');
subplot(235), imshow(a4);
title('(Dylacja - Erozja)/2');

figure;
subplot(231), imshow(b);
title('Peppers');
subplot(232), imshow(b1);
title('Dylacja - Obraz');
subplot(233), imshow(b2);
title('Obraz - Erozja');
subplot(234), imshow(b3);
title('Dylacja - Erozja');
subplot(235), imshow(b4);
title('(Dylacja - Erozja)/2');

figure;
subplot(231), imshow(c);
title('Circles');
subplot(232), imshow(c1);
title('Dylacja - Obraz');
subplot(233), imshow(c2);
title('Obraz - Erozja');
subplot(234), imshow(c3);
title('Dylacja - Erozja');
subplot(235), imshow(c4);
title('(Dylacja - Erozja)/2');

%% Laboratoria 4 - rekonstrukcja morfologiczna
% rekonstrukcja jest cykliczna dylacja i nastepnie laczeniem z obrazem
% czyli dylacja & obraz i usuniecie elementow stycznych z brzegiem
% 1. marker - obraz logiczny o tym samym rozmiarze co obraz wejsciowy z
% wartosciami false
% 2. skopiowanie skrajnych wierszy i kolumn z obrazu do markera
% 3. dopoki nie bedzie roznicy miedzy iteracjami wykonywac:
% - dylacja markera do markera z SE
% - iloczyn logiczny powiekszonego markera z obrazem wejsciowym
% 4. od obrazu wyjsciowego objecie logicznego wyniku rekonstrukcji
close all; clear; clc;
a = imread('blobs.png');
[Nz, Nx] = size(a);
marker = a;
marker(2:Nz - 1, 2:Nx - 1) = false;
b = a;
SE = ones(3);
while ~isequal(marker, b)
    b = marker;
    marker = imdilate(marker, SE) & a;
end
subplot(131), imshow(a);
marker = a & (~marker);
subplot(132), imshow(marker);
subplot(133), imshow(imclearborder(a));

%% Laboratoria 4 - estymacja odległosci geodezyjnej
% start (w,k) 131, 126 do 172, 179
%1. marker z wartoscia true w starcie
%2. rekonstrukcja dopoki w stopie jest false
%3. zliczamy ilosci iteracji koniecznych do zmiany F->T
close all; clear; clc;
a = imread('circles.png');
SE1 = ones(3);
SE2 = [0 1 0; 1 1 1; 0 1 0];
[Nz, Nx] = size(a);
marker1 = false(Nz, Nx);
marker1(131, 126) = true;
marker2 = marker1;
N1 = 0;
while ~marker1(172, 179)
    marker1 = imdilate(marker1, SE1) & a;
    N1 = N1 + 1;
end
N2 = 0;
while ~marker2(172, 179)
    marker2 = imdilate(marker2, SE2) & a;
    N2 = N2 + 1;
end
disp(N1)
disp(N2)
subplot(131), imshow(a);
subplot(132), imshow(marker1);
subplot(133), imshow(marker2);

%% Laboratoria 4 - wypełnianie dziur w obiektach
close all; clear; clc;
a = imread('circles.png');
b = ~a;
c = imclearborder(b);
d = c | a;
subplot(221), imshow(a);
subplot(222), imshow(b), title('Negatyw');
subplot(223), imshow(c), title('Negatyw i wyczysczone brzegi');
subplot(224), imshow(d), title('Bez dziur');

e = imfill(a, 'holes');
figure;
subplot(121), imshow(a);
subplot(122), imshow(e);
% jest tez imreconstruct ale tego nie używaliśmy

%% Laboratoria 5 - hit or miss
% sa rozne oznaczenia:
% - 1 musi byc true (w MatLabie 1)
% - x nie ma znaczenia (w MatLabie 0)
% - 0 musi byc false (w Matlabie -1)
% robimy obraz czarny z odwroconym T w srodku 200x200
close all; clear; clc;
a = false(200, 200);
a(131:160, 41:160) = true;
a(41:130, 86:115) = true;
SE1 = [1 1 -1; 1 1 1; 1 1 1];
hom1 = bwhitmiss(a, SE1) | bwhitmiss(a, rot90(SE1));
subplot(121), imshow(a);
subplot(122), imshow(hom1);

%% Laboratoria 5 - hit or miss do tworzenia figury wypukłej na obiekcie
close all; clear; clc;
a = false(200, 200);
a(131:160, 41:160) = true;
a(41:130, 86:115) = true;
SE1 = [1 1 0; 1 -1 0; 1 0 -1];
SE2 = [1 1 1; 1 -1 0; 0 -1 0];
b = false(200);
subplot(121), imshow(a);
while ~isequal(a, b)
    b = a;
    for k = 1:4
       a = a | bwhitmiss(a, SE1);
       a = bwmorph(a, 'clean');
       a = a | bwhitmiss(a, SE2);
       SE1 = rot90(SE1);
       SE2 = rot90(SE2);
    end
end
subplot(122), imshow(a);

%% Laboratoria 5 - rózne argumenty bwmorph
% bwmorph moze przyjąć argumenty:
% - clean -> fill
% - dilate
% - erode
% - bridge
% - hbreak
% - thin
% - thicken
% - skel 
% - spur
% szkieletowanie jest to zbior wszystkich srodków kół, które
% a) sa w calosci wewnatrz figury
% b) maja minimum dwa punkty wspolne z brzegiem
% zastosowanie ocr:
% - znajdowanie kształtów
% - znajdowanie drogi w labiryncie
% - rozpoznawanie tekstu
close all; clear; clc;
a = imread('circles.png');
b = bwmorph(a, 'clean', Inf);
c = bwmorph(a, 'dilate', Inf);
d = bwmorph(a, 'erode', Inf);
e = bwmorph(a, 'bridge', Inf);
f = bwmorph(a, 'hbreak', Inf);
g = bwmorph(a, 'thin', Inf);
h = bwmorph(a, 'thicken', Inf);
i = bwmorph(a, 'skel', Inf);
j = bwmorph(a, 'spur', Inf);
subplot(341), imshow(a);
title('Circles');
subplot(342), imshow(b);
title('Clean');
subplot(343), imshow(c);
title('Dilate');
subplot(344), imshow(d);
title('Erode');
subplot(345), imshow(e);
title('Bridge');
subplot(346), imshow(f);
title('Hbreak');
subplot(347), imshow(g);
title('Thin');
subplot(348), imshow(h);
title('Thicken');
subplot(349), imshow(i);
title('Skel');
subplot(3,4,10), imshow(j);
title('Spur');

%% Laboratoria 5 - Top Hat i Bot Hat
% Top Hat jest to obraz - otwarcie
% - wyroznia jasniejsze obiekty na ciemniejszym tle
% - duzym elementem strukturalnym > 10
% Bot Hat jest to zamkniecie - obraz
% - wyroznia ciemniejsze obiekty na jasniejszym tle
% - duzym elementem strukturalnym > 10

%% Laboratoria 5 - pole koła (matematycznie i na obrazie)
close all; clear; clc;
a = zeros(256);
for kz = 1:256
    for kx = 1:256
        a(kz, kx) = sqrt((kz-128).^2 + (kx-128).^2);
    end
end
subplot(121), imagesc(a), axis image;

%sposob II
x = -127:128;
[KZ, KX] = meshgrid(x, x);
b = sqrt(KZ.^2 + KX.^2);
subplot(122), imagesc(b), axis image;

R = 5:5:100;
pole_mat = pi*R.^2;
pole_obraz = zeros(size(R));
for k = 1:20
    tt = (a <= R(k));
    pole_obraz(k) = sum(tt(:));
end
figure;
plot(R, pole_mat, '*k', R, pole_obraz, 'r');

%% Laboratoria 5 - obwod koła (matematycznie i na obrazie (filtr górnoprzepustowy i gradient morfologiczny))
close all; clear; clc
a = zeros(256);
for kz = 1:256
    for kx = 1:256
        a(kz, kx) = sqrt((kz-128)^2 + (kx-128)^2);
    end
end
R = 5:5:100;
obw_mat = pi*2*R;
obw_grad = zeros(size(R));
obw_filtr = zeros(size(R));
obw_bwarea = zeros(size(R));
for k = 1:20
    tt = (a <= R(k));
    grad = imdilate(tt, ones(3)) - imerode(tt, ones(3));
    obw_grad(k) = sum(grad(:))/2;
    filter = edge(tt, 'canny');
    obw_filtr(k) = sum(filter(:));
    obw_bwarea(k) = bwarea(bwperim(tt));
end
plot(R, obw_mat, '*k', R, obw_grad, 'r', R, obw_filtr, 'g', R, obw_bwarea, 'b');

%% Laboratoria 5 - odległość między punktami
close all; clear; clc;
a = false(100);
for k = 1:3
    x = ceil(100 * rand(1));
    z = ceil(100 * rand(1));
    a(z, x) = true;
end
b = bwdist(a, 'euclidean');      %odleglosc euklidesowa
c = bwdist(a, 'quasi-euclidean'); %odleglosc quasi euklidesowa
d = bwdist(a, 'cityblock');      %norma L1, cityblock, taksowkowa
e = bwdist(a, 'chessboard');     %norma czebyszewa, szachowa
subplot(231), imshow(a), title('Wylosowane punkty');
subplot(232), imagesc(b), axis image, colorbar('vertical'), title('Euklidesowa');
subplot(233), imagesc(c), axis image, colorbar('vertical'), title('Quasi Euklidesowa');
subplot(234), imagesc(d), axis image, colorbar('vertical'), title('CityBlock');
subplot(235), imagesc(e), axis image, colorbar('vertical'), title('Chessboard');

%% Laboratoria 5 - znajdowanie terenu
% na obrazie new_map.bmp, znalezc miejsce ktore:
% - jest w lesie
% - odległość od drogi głownej jest większa niż 15px
% - odległość od drogi pobocznej jest mniejsza niż 10px
% - odleglość od wody jest większa lub równa 20px
close all; clear; clc;
a = imread('./Obrazy/new_map.bmp');
[Nz, Nx, n] = size(a);
las = a(:,:,1) == 185 & a(:,:,2) == 215 & a(:,:,3) == 170;
droga_glowna = a(:,:,1) == 255 & a(:,:,2) == 245 & a(:,:,3) == 120;
droga_glowna = imclose(droga_glowna, ones(1, 3));
droga_poboczna = a(:,:,1) == 255 & a(:,:,2) == 255 & a(:,:,3) == 255;
droga_poboczna = imopen(droga_poboczna, ones(3));
droga_poboczna = imclose(droga_poboczna, ones(1,3));
woda = a(:,:,1) >= 60 & a(:,:,1) <= 160 & a(:,:,2) >= 150 & a(:,:,2) <= 200 & a(:,:,3) >= 180 & a(:,:,3) <= 255;
wynik = las & (bwdist(droga_glowna) > 15) & (bwdist(droga_poboczna) <= 10) & (bwdist(woda) >= 20);
wynik = imoverlay(a, wynik, 'r');
subplot(231), imshow(a), title('New Map');
subplot(232), imshow(las), title("Teren Lasu");
subplot(233), imshow(droga_glowna), title("Droga Głowna");
subplot(234), imshow(droga_poboczna), title("Droga Poboczna");
subplot(235), imshow(woda), title('Woda');
subplot(236), imshow(wynik), title("Wynik");

%% Laboratoria 6 - były w c++, więc do projektu pomocne, a do kolosa nie

%% Laboratoria 7 - watershed wewnątrz obiektu
close all; clear; clc;
a = false(100, 200);
a(50, [72, 128]) = true;
a = (bwdist(a) <= 30);
D = -bwdist(~a);
L = watershed(D);
b = a & (L > 0);
subplot(221), imshow(a);
subplot(222), imagesc(D), axis image;
subplot(223), imagesc(L), axis image;
subplot(224), imshow(b);

%% Laboratoria 7 - watershed na zewnątrz
close all; clear; clc;
a = false(100, 200);
a(50, [72, 128]) = true;
a = (bwdist(a) <= 30);
temp = imerode(a, ones(23));
D = bwdist(temp);
L = watershed(D);
b = a & (L > 0);
subplot(221), imshow(temp);
subplot(222), imagesc(D), axis image;
subplot(223), imagesc(L), axis image;
subplot(224), imshow(b);

%% Laboratoria 7 - etykietowanie (liczenia pól i obwodów monet)
% wczytac obraz coins.png i zbinearyzowac
close all; clear; clc;
a = imread('coins.png');
bin = a > 90;
bin = medfilt2(bin, [3 3]);
subplot(121), imshow(a);
subplot(122), imshow(bin);
[aseg, N] = bwlabel(bin);
figure;
imagesc(aseg), axis image;
pole = zeros(N, 1);
obwod = zeros(N, 1);
for k = 1:N
    temp = (aseg == k);
    pole(k) = sum(temp(:));
    obwod(k) = bwarea(bwperim(temp));
end
disp(pole)
disp(obwod)

%% Laboratoria 7 - szukanie monet po polu
% na podstawie coins.png, stworzyc nowy obraz, który
% zawiera 5 największych monet w kolorze naturalnym, a pozostałe białe
close all; clear; clc;
a = imread('coins.png');
bin = a > 90;
bin = medfilt2(bin, [3 3]);
subplot(121), imshow(a);
subplot(122), imshow(bin);
[aseg, N] = bwlabel(bin);

figure;
imagesc(aseg), axis image;
pole = zeros(N, 1);
obwod = zeros(N, 1);
for k = 1:N
    temp = (aseg == k);
    pole(k) = sum(temp(:));
    obwod(k) = bwarea(bwperim(temp));
end
prog = median(pole);
wynik = zeros(size(a), 'uint8');
for k = 1:N
    if pole(k, 1) > prog
        wynik = wynik + uint8(aseg == k).*a;
    else
        wynik = wynik + 255 * uint8(aseg == k);
    end
end

figure;
imshow(wynik);

%% Laboratoria 7 - analiza obrazu (schemat i przykład analizy)
% 1. Akwizycja
% 2. Przetwarzanie wstępne:
% - filtracje
% - korekty geometryczne
% - odszumianie
% - przycinanie
% 3. Segmentacja
% 4. Analiza -> przypisanie wartości liczbowych
% 5. Wizualizacja
% Zadanie
% - policzyć pole i obwód każdego z ziaren
% - wyświetlenie histogramów pól i obwodów
close all; clear; clc;
a = imread('rice.png');
b = imtophat(a, strel('disk', 10));
c = imadjust(b);
d = c > 95;
e = bwareaopen(d, 10, 4);
f = imopen(d, ones(5));
temp = imerode(e, ones(5));
D = bwdist(temp);
L = watershed(D);
g = e & (L > 0);
h = imclearborder(g);

subplot(331), imshow(a), title('Rice');
subplot(332), imshow(b), title('Top Hat');
subplot(333), imshow(c), title('Adjust');
subplot(334), imshow(d), title('Binearyzacja');
subplot(335), imshow(e), title('Otwarcie - bwareaopen');
subplot(336), imshow(f), title('Otwarcie - imopen');
subplot(337), imshow(g), title("Watershed");
subplot(338), imshow(h), title("Clear border");

jk = uint8(h).*a;
kj = uint8(~h).*a;
figure;
subplot(121), imshow(jk);
subplot(122), imshow(kj);

[aseg, N] = bwlabel(h);
pole = zeros(N, 1);
obwod = zeros(N, 1);
for k = 1:N
    temp = aseg == k;
    pole(k, 1) = sum(temp(:));
    obwod(k, 1) = bwarea(bwperim(temp));
end
figure;
subplot(121), hist(pole), title("Pola");
subplot(122), hist(obwod), title('Obwody');

%% Laboratoria 8 - transformacja Fouriera
% wykonuje się ją na doublach
close all; clear; clc;
a = imread('cameraman.tif');
a = double(a) / 255;
A = fftshift(fft2(a));
WA = abs(A);
[Nz, Nx] = size(a);
fz = linspace(-0.5, 0.5, Nz);    %*Nz 
fx = linspace(-0.5, 0.5, Nx);    %*Nx
subplot(121), imshow(a);
subplot(122), imagesc(fz, fx, log(WA) * 0.01), axis image;


%% Laboratoria 8 - filtry idealne LP i odwrotna filtracja
close all; clear; clc;
a = imread('cameraman.tif');
a = double(a) / 255;
A = fftshift(fft2(a));
WA = abs(a);
f0 = [0.05, 0.1, 0.2, 0.5];
[Nz, Nx] = size(a);
fx = linspace(-0.5, 0.5, Nx);
fz = linspace(-0.5, 0.5, Nz);
[FX, FZ] = meshgrid(fx, fz);
f = sqrt(FX.^2 + FZ.^2); %odległość od środka
for k = 1 : 4
    LP = f <= f0(k);
    an = real(ifft2(ifftshift(A.*LP)));
    subplot(2, 2, k), imshow(an);
end

%% Laboratoria 8 - filtracja ciągłą nieidealna dolnoprzepustowa, filtracaj Butterwotha
%im wieksze N tym filtr bardziej stromy, 
%tworzymy buttewortha dla N = 2
%jest słabsze ringowanie niż dla idealnego, ale większe rozmycie
close all; clear; clc;
a = imread('cameraman.tif');
a = double(a) / 255;
A = fftshift(fft2(a));
WA = abs(A);
[Nz, Nx] = size(a);
fx = linspace(-0.5, 0.5, Nx);
fz = linspace(-0.5, 0.5, Nz);
[FZ, FX] = meshgrid(fz, fx);
f = sqrt(FZ.^2 + FX.^2);
f0 = [0.05, 0.1, 0.2, 0.5];
for k = 1 : 4
    BW = 1 ./ (1 + (f/f0(k)) .^ 4);
    an = real(ifft2(ifftshift(A .* BW)));
    subplot(2, 2, k), imshow(an);
end

%% Laboratoria 8 - filtracja wycinanie pasma
%dla kazdej palety osobno, R, G, B 

close all; clear; clc
a = imread('./Obrazy/F_dzieciol.png');
a = double(a)/255;
r = a(:,:,1);
g = a(:,:,2);
b = a(:,:,3);
AR = fftshift(fft2(r));
AG = fftshift(fft2(g));
AB = fftshift(fft2(b));
WAR = abs(AR);
WAG = abs(AG);
WAB = abs(AB);
[Nz, Nx, k] = size(a);
fx = linspace(-0.5, 0.5, Nx);
fz = linspace(-0.5, 0.5, Nz);
subplot(221), imshow(a), title('Obraz dzieciol');
subplot(222), imagesc(fx, fz, log(WAR + 0.01)), title('Widmo amplitudowe czerwony');
subplot(223), imagesc(fx, fz, log(WAG + 0.01)), title('Widmo amplitudowe zielona');
subplot(224), imagesc(fx, fz, log(WAB + 0.01)), title('Widmo amplitudowe niebieska');
[FX, FZ] = meshgrid(fx, fz);
BS = abs(FX)>0.17 & abs(FX)<0.24 & abs(FZ)>0.13 & abs(FZ)<0.24;
BS = 1 - BS;
b = a;
figure;
for k = 1 : 3
    A = fftshift(fft2(a(:,:,k)));
    WA = abs(A) .* BS;
    subplot(2,2,k+1), imagesc(fx, fz, log(WA+0.01));
    b(:,:,k) = real(ifft2(ifftshift(BS .* A)));
end
figure;
imshow(b);

%% Laboratoria 8 - korelacja w domenie częstotliwości, szukanie liter 'a'
close all; clear; clc;
bw = imread('text.png');
a = bw(32:45, 88:98);
C = real(ifft2(fft2(bw) .* fft2(rot90(a, 2), 256, 256)));
cmax = max(C(:));
wynik_1 = C > 0.95*cmax;
wynik_2 = imdilate(wynik_1, ones(5, 3));
wynik_3 = imreconstruct(wynik_2, bw);
subplot(221), imshow(bw), title('Tekst');
subplot(222), imshow(wynik_1), title('Znalezione "a", przed dylacja');
subplot(223), imshow(wynik_2), title('Znalezione "a", po dylacji');
subplot(224), imshow(wynik_3), title('Znalezione "a", po rekonstrukcji');

%% Laboratoria 8 - korelacja w domenie częstotliwości, szukanie liter 'r'
close all; clear; clc;
bw = imread('text.png');
a = bw(32:45, 4:13);
C1 = real(ifft2(fft2(bw) .* fft2(rot90(a, 2), 256, 256)));
C2 = real(ifft2(fft2(~bw) .* fft2(rot90(~a, 2), 256, 256)));
C = C1 + C2;
cmax = max(C(:));
wynik_1 = C > 0.99 * cmax;
wynik_2 = circshift(wynik_1, [0, -4]);
wynik_3 = imdilate(wynik_2, ones(5, 3));
wynik_4 = imreconstruct(wynik_3, bw);
subplot(221), imshow(wynik_1), title('Znalezione "r"');
subplot(222), imshow(wynik_2), title('Znalezione "r", po circshift')
subplot(223), imshow(wynik_3), title('Znalezione "r", po dylacji');
subplot(224), imshow(wynik_4), title('Znalezione "r", po rekonstrukcji');

%% Laboratoria 8 - transformacja falkowa, czasowo częstotliwościwa
%wavemenu w konsoli
close all; clear; clc;
a = imread('cameraman.tif');
[C, L] = wavedec2(a, 2, 'sym3');
L2 = L(:,1) .* L(:,2);
A2 = C(1:L2(1));
H2 = C(L2(1) + 1:L2(1) + L2(2));
V2 = C(2*L2(1) + 1:L2(1) + 2*L2(2));
D2 = C(3*L2(1) + 1:L2(1) + 3*L2(2));
H1 = C(4*L2(1) + 1:4*L2(1) + L2(3));
V1 = C(4*L2(1) + L2(3) + 1:4*L2(1) + 2*L2(3));
D1 = C(4*L2(1) + 2*L2(3) + 1:4*L2(1) + 3*L2(3));
temp = reshape(D1, L(3,:));
imagesc(temp), axis image;

%% Laboratoria 8 - transformacja falkowa, usuwanie składowych, zadanko
%korzystajac z rekonstrukcji wyzerowac A2, H1, D1, waverec2(C1, L, 'sym3')
% i porownac z obrazem wejsciowym
close all; clear; clc;
a = imread('cameraman.tif');
[C, L] = wavedec2(a, 2, 'sym3');
C1 = C;
C2 = C;
C3 = C;
L2 = L(:, 1).*L(:,2);
C1(1:L2(1)) = 0;
an_1 = waverec2(C1, L, 'sym3');

C2(4*L2(1) + 1:4*L2(1) + L2(3)) = 0;
an_2 = waverec2(C2, L, 'sym3');

C3(4*L2(1) + 2*L2(3) + 1:4*L2(1) + 3*L2(3)) = 0;
an_3 = waverec2(C3, L, 'sym3');
subplot(221), imagesc(a), axis image, title('Cameraman');
subplot(222), imagesc(an_1), axis image, title('Rekonstrukcja z wyzerowanym A2');
subplot(223), imagesc(an_2), axis image, title('Rekonstrukcja z wyzerowanym H1');
subplot(224), imagesc(an_3), axis image, title('Rekonstrukcja z wyzerowanym D1');

%% Laboratoria 9 - transformata Gabora
% filtry zespolone
% iloczyn kartezjanski dlugosci fali, ilosci katow, orientacji funkcji
% gaussa, odlgelosci od srodka
close all; clear; clc;
a = imread('./Obrazy/Gabor.png');
a = imresize(a, 0.5); % zmniejszenie obrazu dwukrotnie
dlug_fali = 2.^(1:6); % okres sinusoidy (1:6)
krok = 22.5; % o jaki stopien rotowac
katy = 0 : krok : 180 - krok; % kolejne rotacje o krok
g = gabor(dlug_fali, katy); % bank
magn = imgaborfilt(a, g); % magnituda, konwolucja obrazu z filtrami
subplot(131), imshow(a), title('Obraz wejściowy');
subplot(132), imagesc(abs(g(1).SpatialKernel)), title('Spatial Kernel g(1)');
subplot(133), imagesc(abs(g(48).SpatialKernel)), title('Spatial Kernel g(48)');

figure;
subplot(121), imagesc(magn(:, :, 1)), title('magnituda 1');
subplot(122), imagesc(magn(:, :, 47)), title('magnituda 47');

for k = 1 : length(g)
   odch = 1.5 * g(k).Wavelength;
   magn(:, :, k) = imgaussfilt(magn(:, :, k), odch);
end

%% Laboratoria 9 - kmeans podział na dwa klastry
[Nz, Nx] = size(a);
x = 1 : Nx;
z = 1 : Nz;
[XX, ZZ] = meshgrid(x, z);
zbior = cat(3, magn, XX);
zbior = cat(3, zbior, ZZ);
D2 = reshape(zbior, Nx * Nz, []);
D1 = D2 - mean(D2);
D1 = D1 ./ std(D1);
L = kmeans(D1, 2, 'Replicates', 5);
wynik = reshape(L, [Nz, Nx]);
imagesc(wynik);

%% Laboratoria 9 - działanie dalej na klastrze z tekstem, usuniecie dziur, czyszczenie brzegu, usuniecie obiektow o polu mniejszym niz 600 pikseli
bw = (wynik == 1);
bw1 = imclearborder(bw);
bw2 = bwareaopen(bw1, 600);
bw3 = imfill(bw2, 'holes');
an = imoverlay(a, ~bw3, 'k');
subplot(231), imagesc(wynik), title('Przed działaniem');
subplot(232), imshow(bw), title('Binearyzacja');
subplot(233), imshow(bw1), title('Czyszczenie brzegow');
subplot(234), imshow(bw2), title('Usuniecie obiektow o polu < 600');
subplot(235), imshow(bw3), title('Usuniecie dzior');
subplot(236), imshow(an), title('Koncowy wynik');

%% Laboratoria 9 - zadanie
% zmieniac dlugosc fali oraz krok, stworzyc taka segmentacje aby kazdy
% wyraz byl osobno a nie razem jak teraz
% trzeba bylo zmniejszyc dlugosc fali i zagescic kroki

%% Laboratoria 9 - transformata cosinusowa
% o rozmiarze obrazka, lewy gorny rog transformaty 
% funkcja dct2, oraz idct2 dla odwrotnej
% magnituda zamiast WA i WF poniewaz jestesmy w liczbach rzeczywistych
close all; clear; clc;
a = imread('cameraman.tif');
a = double(a) / 255;
A = dct2(a);
subplot(121), imshow(a), title("Cameraman");
subplot(122), imagesc(log(abs(A) + 0.01)) , title('Transformata cosinusowa');

%odwrotna transformata idealna
r = [5, 10, 25, 50, 100];
odleglosc = zeros(size(a));
odleglosc(1,1) = 1;
odleglosc = bwdist(odleglosc);
figure;
subplot(231), imagesc(odleglosc), title('Mapa odległości');
for k = 1:length(r)
    maska = odleglosc <= r(k);
    an = idct2(A.*maska);
    str = ['R = ', num2str(r(k))];
    subplot(2, 3, k + 1), imshow(an), title(str)
end

figure;
th = [0.01 0.05 0.1 0.5 1.0];
for k = 1 : length(th)
    maska2 = abs(A) >= th(k);
    an2 = idct2(A.*maska2);
    disp(sum(maska2(:) == 0));
    subplot(2, 3, k), imshow(an2), title(['th = ', num2str(th(k))]);
end
subplot(236), imshow(a);

%% Laboratoria 9 - konwersja JPEG92
%kompresja i dekompresja
% 1. zmiana palety barw z RGB na YCbCr
% 2. przesuwamy srednia do zera
% 3. downsampling, usuwamy co 2-3 wiersz z palet Cb, Cr
% 4. dzielimy warstwy na bloki 8x8 i dla kazdego bloku
% a) transformate cosinusowa
% b) podział przez tablice kwantyfikujące Qy i Qc
% c) zaokraglenie
% b i c daja razem kwantyfikacje
% 5. zigzak - zamiana bloku 8 na 8 na wektor w specyficzny sposob, w taki
% 6. kompresja znana z winrara itp
% 7. zapis do pliku
% sposob ze szybko dochodzimy do ostatniej niezerowej wartosci
close all; clear; clc;
Qy = [16 11 10 16 24 40 51 61
12 12 14 19 26 58 60 55
14 13 16 24 40 57 69 56
14 17 22 29 51 87 80 62
18 22 37 56 68 109 103 77
24 35 55 64 81 104 113 92
49 64 78 87 103 121 120 101
72 92 95 98 112 100 103 99];

Qc = [17 18 24 47 99 99 99 99
18 21 26 66 99 99 99 99
24 26 56 99 99 99 99 99
47 69 99 99 99 99 99 99
99 99 99 99 99 99 99 99
99 99 99 99 99 99 99 99
99 99 99 99 99 99 99 99
99 99 99 99 99 99 99 99];
skala = 1.1;
Qc = skala * Qc;
Qy = skala * Qy;
a = imread('peppers.png');
b = rgb2ycbcr(a);
b = double(b);
b = b - 128;  % mean(b(:)) tylko wtedy trzeba zapisac do zmiennej
[Nz, Nx, k] = size(a);
ile_zer = 0;
for k = 1 : 3
    for z = 1:8:Nz
        for x = 1:8:Nx
            tt = b(z:z + 7, x:x + 7, k);
            tt = dct2(tt);
            if k == 1
                tt = tt ./ Qy;
            else
                tt = tt ./ Qc;
            end
            tt = round(tt);
            ile_zer = ile_zer + sum(tt(:) == 0);
            if k == 1
                tt = tt .* Qy;
            else
                tt = tt .* Qc;
            end
            tt = idct2(tt);
            b(z: z + 7, x:x + 7, k) = tt;
        end
    end
end
b = uint8(b + 128);
b = ycbcr2rgb(b);
disp(100*ile_zer/(Nz * Nx * k))
subplot(121), imshow(a);
subplot(122), imshow(b);

%% Laboratoria 9 - transformata Radona
% transformata nie jest odwrotna
close all; clear; clc;
a = false(120, 200);
a(31:90, 51:150) = true;
subplot(131), imshow(a), title('Biały prostokąt na czarnym tle');
kat = 0:1:180;
[R, X] = radon(a, kat);
subplot(132), imagesc(kat, X, R), title("Transformata Radona"), colormap(hot(256));
an = iradon(R, kat);
subplot(133), imshow(an), title("Odwrotona do transformaty Radona");

%% Laboratoria 10 - transformata Hougha, znajdowanie linii
close all; clear; clc;
a = imread('blobs.png');
[H, T, R] = hough(a); %H - obraz transformaty, T - theta, R - rho (czyli osie)
subplot(121), imshow(a), title('Blobs.png');
subplot(122), imagesc(T, R, H), title('Transformata Hougha');
pik = houghpeaks(H, 10); %ile linii chcemy odzyskać jako drugi argument
L = houghlines(a, T, R, pik, 'FillGap', 5);

figure;
imshow(a), title('Blobs.png'); hold on;
kk = 1;
max = sqrt(abs(L(1).point2(1) - L(1).point1(1)) ^ 2 + abs(L(1).point2(2) - L(1).point1(2)) ^ 2);
for k = 1 : 10
    temp = sqrt(abs(L(k).point2(1) - L(k).point1(1)) ^ 2 + abs(L(k).point2(2) - L(k).point1(2)) ^ 2);
    if temp > max
        kk = k;
    end
    line([L(k).point1(1), L(k).point2(1)], [L(k).point1(2), L(k).point2(2)], 'color', 'red');
end
line([L(kk).point1(1), L(kk).point2(1)], [L(kk).point1(2), L(kk).point2(2)], 'color', 'green');
hold off;

%% Laboratoria 10 - Analiza Obrazu
close all; clear; clc;
[map, leg] = imread('./Obrazy/w_shape.png'); %mapa i legenda - obraz indeksowany
a = ind2rgb(map, leg);
imagesc(map), title('Sprawdzanie kolorow');
bin = map ~= 11;
figure;
imshow(bin), title('Binearyzacja');
[aseg, N] = bwlabel(bin);
rp = regionprops(aseg, 'all'); %sprawdzamy rózne statystyki obiektów

%% Laboratoria 10 - Analiza Obrazów, szukanie obiektów
kolo = false(size(bin));
kwadrat = false(size(bin));
elipsa = false(size(bin));
gwiazdy = false(size(bin));
BWK = 0;
for k = 1 : N
    temp = aseg == k;
    obwod = bwarea(bwperim(temp));
    pole = bwarea(temp);
    bwk = 4 * pi * pole / (obwod * obwod);
    BWK = BWK + bwk * temp;
    if abs(bwk - 1) < 0.05
        kolo = kolo | temp;
    end
    if abs(bwk - pi/4) < 0.02
        kwadrat = kwadrat | temp;
    end
    pole_mat = pi * rp(k).MinorAxisLength * rp(k).MajorAxisLength/4;
    if(abs(pole_mat/pole - 1)) < 0.02
        elipsa = elipsa | temp;
    end
    con = rp(k).ConvexArea - rp(k).Area;
    if bwk <= 0.27 & bwk >= 0.24 & rp(k).EulerNumber == 1  
        gwiazdy = gwiazdy | temp;
    end
end
figure;
subplot(231), imshow(a), title('Wejsciowy obraz');
subplot(232), imshow(kolo), title('Znalezione kola');
subplot(233), imshow(kwadrat), title('Znalezione kwadraty');
subplot(234), imshow(elipsa), title('Znalezione elipsy');
subplot(235), imshow(gwiazdy), title('Znalezione gwiazdy');
subplot(236), imagesc(BWK), title("Rozklad BWK");

%% Laboratoria 10 - analiza wykresu
close all; clear; clc;
a = imread('./Obrazy/wykres.png');
subplot(121), imshow(a), title('Wykres');
y_min = 113;
y_max = 1312;
x_min = 346;
x_max = 2375;
bin = a(y_min:y_max, x_min:x_max, 1) == 126;
subplot(122), imshow(bin), title('Interesujaca nas czesc wykresu');
[Nz, Nx] = size(bin);
t = zeros(Nx, 1);
p = zeros(Nx, 1);
k = 1;
for kx = 1 : Nx
    suma = sum(bin(:, kx));
    if suma > 0
        k = k + 1;
        pmin = find(bin(:, kx), 1, 'first');
        pmax = find(bin(:, kx), 1, 'last');
        p(k) = pmax - pmin / 2;
        t(k) = kx;
    end
end
t(k : end) = [];
p(k : end) = [];
t = 350 * t / Nx;
p = 35000 - 35000 * p / Nz;
t2 = 0 : 2 : 24 * 14;
p2 = interp1(t, p, t2);
figure;
plot(t, p, 'r', t2, p2, 'go')