% Rodzial 17 Zielniski
%                       Mateusz Krupnik
% �wiczenie: Wykorzystanie transformacji 2D DFT i 2D DCT do filtracji obrazu
clear all; close all; clc;
% Inicjalizacja ? wczytaj obraz
[x,cmap] = imread('cameraman.tif');     % wczytaj obraz do "x" i jego palet� kolor�w do "cmap"
imshow(x,cmap), title('Obraz');	% poka� obraz wykorzystuj�c jego palet�
[M, N] = size(x);       % odczytaj liczb� wierszy i kolumn; za�o�enie M=N !!!
x = im2double(x);           % zamie� reprezentacj� pikseli

% Macierz transformacji 1D DCT oraz 1D DFT
n=0:N-1; % numery pr�bek funkcji bazowych
c = [sqrt(1/N) sqrt(2/N)*ones(1,N-1)]; f = 1/sqrt(N); % wsp�czynniki normalizuj�ce
for k=0:N-1                 % wyznacz macierz transformacji
    C(k+1,n+1) = c(k+1) * cos( pi*k*(n+1/2) / N ); % funkcje bazowe 1D DCT
    F(k+1,n+1) = f * exp(j*2*pi/N*k*n); % funkcje bazowe 1D DFT
end

% JEDNA LINIA - transformacja DFT i DCT wiersza
Nr = 100; K = 2; y = x;	% numer linii, szeroko�� znacznika, kopia obrazu
linia = x(Nr,1:N);      % pobranie linii
y(Nr-K:Nr+K,1:N) = 0*ones(2*K+1,N);	% zaznacz wybran� lini� czarnym kolorem
% poka� obraz z czarn� lini�
figure(2)
subplot(221); imshow(y,cmap); title('Obraz');
% poka� wykres linii obrazu
subplot(222); plot(linia); title('Jedna linia'); 
% DFT linii obrazu
subplot(223); plot( abs(fft(linia))/N ); title('|DFT|');
grid;
% DCT (Matlab) linii
% subplot(212); plot( dct(linia)/sqrt(N) ); title('DCT'); 
 % DCT (nasze) linii
subplot(224); plot( (conj(C)*linia')/sqrt(N) ); title('DCT');
grid;

% FILTRACJA za pomoc� 2D DCT
% maska cz�stotliwo�ciowa, lewy g�rny r�g
K = 64; H = zeros(M,N); H(1:K,1:K) = ones(K,K);
% X = dct2(x); % 2D DCT ? MATLABA
X = conj(C) * x * conj(C).'; % 2D DCT ? NASZE (dla DCT conj(C)=C)
Y = X .* H; % iloczyn widma DCT i maski
% y = idct2(Y); % 2D IDCT ? MATLABA
y = C.' * Y * C; % 2D IDCT ? NASZE
XdB = skaladB( X ); YdB = skaladB( Y ); % wyskalowanie intensywno�ci pikseli w dB
figure(3)
subplot(221); imshow(x, cmap); title('Obraz'); % dalej tylko wizualizacja
subplot(222); imshow(XdB, cmap); title('Widmo DCT');
subplot(223); imshow(YdB, cmap); title('Widmo DCT + Maska');
subplot(224); imshow(y(1:128,65:192), cmap); title('Fragment wyniku');

% FILTRACJA 2D za pomoc� 2D DFT (fftshift2D - przestawianie �wiartek widma 2D DFT, patrz rys. 22.11a)
% maska cz�stotliwo�ciowa H (MxN)
K = 32; H = zeros(M,N);
% �rodek = (M/2+1, N/2+1)
H(M/2+1-K : M/2+1+K, N/2+1-K : N/2+1+K) = ones(2*K+1,2*K+1);
h = fftshift2D( real( ifft2( fftshift2D(H) ) ) ); % odpowied� impulsowa maski
figure(4)
subplot(121); imshow(255*H,cmap); title('Maska Freq'); % rysunek maski
subplot(122); mesh(h); title('Odpowied� impulsowa'); % rysunek jej odp. impulsowej


% X = fft2(x)/N; % transformacja Fouriera 2D DFT ? MATLABA
X = conj(F) * x * conj(F).';	% transformacja Fouriera 2D DFT ? NASZA
Xp = fftshift2D(X);             % przestawienie miejscami �wiartek widma
Yp = Xp .* H;               % filtracja = iloczyn widma 2D DFT i maski 2D
Y = fftshift2D(Yp);         % powrotne przestawienie �wiartek widma
% y1 = ifft2(Y)*N;          % odwrotna transformacja Fouriera 2D IDFT ? MATLABA
y1 = F.' * Y * F;           % odwrotna transformacja Fouriera 2D IDFT ? NASZA
y1 = real(y1);              % cz�� rzeczywista, urojona r�wna zeru
y1f = y1(1:128,65:192);     % wybranie fragmentu obrazu do wizualizacji

XdB = skaladB( X ); XpdB = skaladB( Xp );
YdB = skaladB( Y ); YpdB = skaladB( Yp );
figure(5)
subplot(231); imshow(x, cmap); title('1. Obraz');
subplot(232); imshow(XdB, cmap); title('2. 2D DFT');
subplot(233); imshow(XpdB, cmap); title('3. Po przestawieniu');
subplot(234); imshow(YpdB, cmap); title('4. Po filtrze');
subplot(235); imshow(YdB, cmap); title('5. Po przestawieniu');
subplot(236); imshow(y1f, cmap); title('6. Fragment wyniku');

% FILTRACJA 2D za pomoc� splotu 2D
L = 32;
y2 = conv2(x, h(M/2+1-L:M/2+1+L, N/2+1-L:N/2+1+L),'same');
figure(6)
subplot(121); imshow(y1,cmap);
title('Obraz po filtrze FREQ - splot cykliczny');
subplot(122); imshow(y2,cmap);
title('Obraz po filtrze CONV - splot liniowy');

%% �wiczenie: Projektowanie filtr�w 2D
clear all;
L = 15; % szeroko�� macierzy wag filtra (nieparzysta: 3, 5, 7, 9, ...)
K = (L-1)/2; df = 0.5/K; % zmienne pomocnicze do generacji wag

m = ones(L,1)*(-K:K); % i opisu osi rysunk�w
n = (-K:K)'*ones(1,L);
fm = ones(L,1)*(-0.5:df:0.5);
fn = (-0.5:df:0.5)'*ones(1,L);

% Wczytaj obraz do filtracji
[x,cmap] = imread('cameraman.tif');
imshow(x,cmap), title('Obraz');
[N, M] = size(x);
x = im2double(x);

% FILTRY POCHODZ�CE OD FUNKCJI GAUSSA
sigma = 1.4; df = 0.5/K;
g0 = 1/(2*pi*sigma^2) .* exp(-(m.^2+n.^2)/(2*sigma^2));	% funkcja Gaussa
g1m = -m/(sigma^2) .* g0;       % pochodna wzgl�dem osi m
g1n = -n/(sigma^2) .* g0;       % pochodna wzgl�dem osi n
g2 = (m.^2 + n.^2 - 2*sigma^2)/(sigma^4) .* g0;	% laplasjan funkcji Gaussa
figure(7); sgtitle('Filtry pochodzace od f. Gaussa');
subplot(221); mesh(m,n,g0); title('Filtr Gaussa');
subplot(222); mesh(m,n,g2); title('Laplasjan f. Gaussa');
colormap([0 0 0]);
subplot(223); imshow( conv2(x,g0,'same'),cmap );
subplot(224); imshow( conv2(x,g2,'same'),cmap );
figure(8); sgtitle('a');
subplot(231); mesh(m,n,g1m); title('Gradient "m"');
subplot(232); mesh(m,n,g1n); title('Gradient "n"');
colormap([0 0 0]);
subplot(234); imshow( conv2(x,g1m,'same'),cmap );
subplot(235); imshow( conv2(x,g1n,'same'),cmap );
subplot(236);
imshow(sqrt(conv2(x,g1m,'same').^2+conv2(x,g1n,'same').^2),cmap);


%% METODA OKIEN
chka = 0; % 0 = charakterystyka prostok�tna, 1 = ko�owa
w = hamming(L); w = w * w'; % okno 2D
figure(9);
subplot(111); mesh(m,n,w); colormap([0 0 0]); title('2D Okno');
for chka=0:1
    if(chka==0) % Charakterystyka prostok�tna - odp. impulsowa dw�ch filtr�w LowPass
        f0=0.25; wc=pi*f0;
        sinc=sin(wc*(-K:K))./(pi.*(-K:K));
        sinc(K+1)=f0; lp1=sinc'*sinc;

        f0=0.50; wc=pi*f0;
        sinc=sin(wc*(-K:K))./(pi.*(-K:K));
        sinc(K+1)=f0; lp2=sinc'*sinc;
        chkat = "prostokatna"
    else % Charakterystyka ko�owa - odp. impulsowa dw�ch filtr�w LowPass
        f0=0.25; wc=pi*f0;
        lp1=wc*besselj( 1,wc*sqrt(m.^2 + n.^2))./(2*pi*sqrt(m.^2+n.^2));
        lp1(K+1,K+1)= wc^2/(4*pi);
        f0=0.50; wc=pi*f0;
        lp2=wc*besselj( 1,wc*sqrt(m.^2+n.^2))./(2*pi*sqrt(m.^2+n.^2) );
        lp2(K+1,K+1)= wc^2/(4*pi);
        chkat = "kolowa"
    end
    lp = lp1; % LowPass bez okna 2D
    lpw = lp .* w; % z oknem
    hp = - lp1; hp(K+1,K+1) = 1 - lp1(K+1,K+1); % HighPass bez okna 2D
    hpw = hp .* w; % z oknem
    bp = lp1 - lp2; % BandPass bez okna 2D
    bpw = bp .* w; % z oknem
    bs = - bp; bs(K+1,K+1) = 1 - bp(K+1,K+1); % BandStop bez okna 2D
    bsw = bs .* w; % z oknem
    for typ = 1 : 4 % poka� odp. impulsow�, jej widmo i przefiltrowany obraz
        switch (typ) % wybierz typ filtra
            case 1, h = lp; hw = lpw; filtr = "LowPass"; % LP
            case 2, h = hp; hw = hpw; filtr = "HighPass"; % HP
            case 3, h = bp; hw = bpw; filtr = "BandPass"; % BP
            case 4, h = bs; hw = bsw; filtr = "BandStop"; % BS
        end
        figure(9+(chka*4)+typ)
        sgtitle("Filtracja: " + chkat + ", " + filtr);
        subplot(321);
        mesh(m,n,h); title('Filtr h(m,n)');
        subplot(322);
        mesh(m,n,hw); title('Filtr hw(m,n)');
        subplot(323);
        mesh(fm,fn,abs( fftshift2D(fft2(h)) ) ); title('|H(fm,fn)|');
        subplot(324);
        mesh(fm,fn,abs( fftshift2D(fft2(hw)) ) ); title('|Hw(fm,fn)|');
        colormap([0 0 0]);
 
        subplot(325); y = conv2(x, h,'same');
        imshow(y,[min(min(y)),max(max(y))]);
        subplot(326); y = conv2(x, hw,'same');
        imshow(y,[min(min(y)),max(max(y))]);
    end
end

%% FILTRY PROJEKTOWANE W DZIEDZINIE CZ�STOTLIWO�CI
% Zmienne pomocnicze do generacji wag i opisu osi rysunk�w
N=L+1; df = 0.5/(K+1);
m = ones(L+1,1)*(-(K+1):K); n = (-(K+1):K)'*ones(1,L+1);
fm = ones(L+1,1)*(-0.5:df:0.5-df); fn = (-0.5:df:0.5-df)'*ones(1,L+1);
% Okno 2D - jego kszta�t i widmo
w = hamming(N); w = w * w';
% Zadana charakterystyka cz�stotliwo�ciowa - ko�owa lub prostok�tna
Q = round(K/2); % szeroko�� filtra LP
H = zeros(N,N);
for k = N/2+1-Q : N/2+1+Q % ko�owa
    for l = N/2+1-Q : N/2+1+Q
        if( (k-N/2-1)^2 + (l-N/2-1)^2 <= Q^2)
            H(k,l) = 1;
        else
            H(k,l) = 0;
        end
    end
end

% H(N/2+1-Q : N/2+1+Q, N/2+1-Q : N/2+1+Q) = ones(2*Q+1,2*Q+1); % prostok�tna
% Zaprojektowanie filtra i sprawdzenie jego dzia�ania
h = real( ifft2(fftshift2D(H)) ); h = fftshift2D(h); % odpowied� impulsowa
hw = h .* w; % odp. impulsowa z oknem
Hw = abs( fftshift2D( fft2( hw ) ) ); % jej widmo
y = conv2(x, hw, 'same'); % filtracja
% Rysunki
figure(18)
sgtitle('Okno filtracji i jego ch. czest.');
subplot(121);
mesh( m,n,w ); title('Okno 2D w(m,n)');
subplot(122);
mesh( fm,fn,abs( fftshift2D(fft2(w)) ) ); title('|W(fm,fn)|');
colormap([0 0 0]);
figure(19);
sgtitle('Filtr projketowany w dz. czest.');
subplot(221); mesh(fm,fn,H); title('Zadane H(fm,fn)');
subplot(222); mesh(m,n,h); title('Filtr 2D h(m,n)');
subplot(223); mesh(fm,fn,Hw); title('Ch-ka |Hw(fm,fn)|');
subplot(224); mesh(m,n,hw); title('Filtr 2D hw(m,n) z oknem');
colormap([0 0 0]);
figure(20);
sgtitle('Filtracja obrazu charaktarystyka projektowana w d. czest.');
subplot(121); imshow(y, cmap);
title('Ca�y obraz po filtrze'); y = y(1:128,65:192);
subplot(122); imshow(y,[min(min(y)),max(max(y))]);
title('Tylko fragment');



%%% FUNKCJE ZAGNIEZDZONE %%%
function Y = fftshift2D( X )
    % przestawianie �wiartek widma 2D DFT
    [M N] = size(X);
    Y(M/2+1:M,N/2+1:N) = X(1:M/2,1:N/2);
    Y(1:M/2,1:N/2) = X(M/2+1:M,N/2+1:N);
    Y(M/2+1:M,1:N/2) = X(1:M/2,N/2+1:N);
    Y(1:M/2,N/2+1:N) = X(M/2+1:M,1:N/2);
end
function XdB = skaladB(X)
    % skalowanie intensywno�ci pikseli obrazu w decybelach
    XdB = log10(abs(X)+1); maxXdB = max(max(XdB)); minXdB = min(min(XdB));
    XdB = (XdB-minXdB)/(maxXdB-minXdB)*255;
end