% Lab 7. Por�wnanie dzia�ania filtra rekursywnego i nierekursywnego. Metody
% fir1 i butter.
%                           Mateusz Krupnik'
clc; clear all; close all;

% Dane sygna��w z lab. 6
f1=100; f2=250; f3=400; fs=1000;    % Czestotliwosci skladowych i Nyquista
A1=1; A2=0.8; A3=0.65;              % Aplitudy skladowych
t=0:(1/fs):1.023;                   % Wektor czasu
% Sygna� wymuszenia - sk�adowa 3 harmonicznych
x=A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t)+A3*sin(2*pi*f3*t);
% Parametry odci�cia
fo1 = 300; fo2 = 200;           % Czestotliwosci ociecia
wn1 = fo1*2/fs; wn2 = fo2*2/fs; % Znormalizowane czestotliwosci

%% Filtr butterwotha dolnoprzepustowy
rzad = 32;          % rz�d filtra
[B, A] = butter(32, wn1);   % dla czestosci odciecia nr 1
x_filtered = filter(B, A, x);   % filtracja syg. przez uk�ad
% Obliczenie widma i mocy za pomoc� funkcji z lab 4.
[f_w, Moc, Wid] = fft_from_signal([x; x_filtered], fs);

% Wykresy
figure(1);
sgtitle(['Filtr Butterworth, rzad=' num2str(rzad) ' f_o=' num2str(fo1)]);
subplot(211);
plot(t, x); title('Sygna� wymuszaj�cy');
xlabel('Czas [s]'); ylabel('Amp.'); grid;
subplot(212);
plot(t, x_filtered); title('Sygna� po filtracji');
xlabel('Czas [s]'); ylabel('Amp.'); grid;

figure(2);
sgtitle(['Filtr Butterworth, rzad=' num2str(rzad) ', f_o=' num2str(fo1)]);
subplot(211);
plot(f_w, Wid(1,:)); title('Widmo sygna�u wymuszaj�cego');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;
subplot(212);
plot(f_w, Wid(2,:)); title('Widmo po filtracji');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;

%% Filtr metoda probkowania w dziedzinie czestotlwiosci
% filtr nierekursywny
rzad = 16;  % rzad filtra
b1 = fir1(rzad, wn1);
x_fir2 = filter(b1, 1, x);
% Obliczenie widma i mocy za pomoc� funkcji z lab 4.
[f_w, Moc, Wid] = fft_from_signal([x; x_fir2], fs);

% Wykresy
figure(3);
sgtitle(['Filtr FIR2, rzad=' num2str(rzad) ' f_o=' num2str(fo1)]);
subplot(211);
plot(t, x); title('Sygna� wymuszaj�cy');
xlabel('Czas [s]'); ylabel('Amp.'); grid;
subplot(212);
plot(t, x_fir2); title('Sygna� po filtracji');
xlabel('Czas [s]'); ylabel('Amp.'); grid;

figure(4);
sgtitle(['Filtr FIR2, rzad=' num2str(rzad) ', f_o=' num2str(fo1)]);
subplot(211);
plot(f_w, Wid(1,:)); title('Widmo sygna�u wymuszaj�cego');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;
subplot(212);
plot(f_w, Wid(2,:)); title('Widmo po filtracji');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;

% Wykresy por�wnawcze filtracji
figure(5);
spectrogram(x, 'yaxis'); title('Sygna� wymuszaj�cy');
figure(6);
spectrogram(x_filtered, 'yaxis'); title('Sygna� po filtracji - IIR');
figure(7);
spectrogram(x_fir2, 'yaxis'); title('Sygna� po filtracji - FIR');

%% Filtrowanie i kr�tkoczasowa analiza cz�stotliwo�ciowa
% sygna��w o zmiennej czestotliwo�ci
T = 0:(1/fs):1.023;
% chirp(wektor czasu, czestestotliwosc start, koniec narastania,
% czestotliwoscdocelowa)
X = chirp(T,50,1.023,450);      % funkcj generuj�ca sygna�
% filtracja przez filtry FIR i IIR
X_fir2 = filter(b1, 1, X);
X_iir = filter(B, A, X);

figure(8)
subplot(311); plot(T, X);
title('Sygna� wymuszaj�cy'); xlabel('Czas [s]'); ylabel('Amp.'); grid;

subplot(312); plot(T, X_fir2);
title('Filtracja FIR'); xlabel('Czas [s]'); ylabel('Amp.'); grid;

subplot(313); plot(T, X_iir);
title('Filtracja IIR'); xlabel('Czas [s]'); ylabel('Amp.'); grid;


figure(9)
% SFFT od sygnalu X, z oknem 256, zak�adk� okien 250, i 256 punktow� FFT
spectrogram(X, 256, 250, 256, 1E3, 'yaxis'); title('Sygna� wymuszaj�cy');

figure(10)
spectrogram(X_fir2, 256, 250, 256, 1E3, 'yaxis');
title('Sygna� po fitlracji FIR');

figure(11)
spectrogram(X_iir, 256, 250, 256, 1E3, 'yaxis');
title('Sygna� po fitlracji IIR');
