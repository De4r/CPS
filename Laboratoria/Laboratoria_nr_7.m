% Lab 7. Porównanie dzia³ania filtra rekursywnego i nierekursywnego. Metody
% fir1 i butter.
%                           Mateusz Krupnik'
clc; clear all; close all;

% Dane sygna³ów z lab. 6
f1=100; f2=250; f3=400; fs=1000;    % Czestotliwosci skladowych i Nyquista
A1=1; A2=0.8; A3=0.65;              % Aplitudy skladowych
t=0:(1/fs):1.023;                   % Wektor czasu
% Sygna³ wymuszenia - sk³adowa 3 harmonicznych
x=A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t)+A3*sin(2*pi*f3*t);
% Parametry odciêcia
fo1 = 300; fo2 = 200;           % Czestotliwosci ociecia
wn1 = fo1*2/fs; wn2 = fo2*2/fs; % Znormalizowane czestotliwosci

%% Filtr dolnoprzepustowy
% Filtr Butterwortha
rzad = 32;          % rz¹d filtra
[B, A] = butter(rzad, wn1);   % dla czestosci odciecia nr 1
x_filtered = filter(B, A, x);   % filtracja syg. przez uk³ad
% Obliczenie widma i mocy za pomoc¹ funkcji z lab 4.
[f_w, Moc, Wid] = fft_from_signal([x; x_filtered], fs);

% Wykresy
figure(1);
sgtitle(['Filtr dolnoprzepustowy Butterworth, rzad=' num2str(rzad)...
    ' f_o=' num2str(fo1)]);
subplot(211);
plot(t, x); title('Sygna³ wymuszaj¹cy');
xlabel('Czas [s]'); ylabel('Amp.'); grid;
subplot(212);
plot(t, x_filtered); title('Sygna³ po filtracji');
xlabel('Czas [s]'); ylabel('Amp.'); grid;

figure(2);
sgtitle(['Filtr dolnoprzepustowy Butterworth, rzad=' num2str(rzad)...
    ', f_o=' num2str(fo1)]);
subplot(211);
plot(f_w, Wid(1,:)); title('Widmo sygna³u wymuszaj¹cego');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;
subplot(212);
plot(f_w, Wid(2,:)); title('Widmo po filtracji');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;

% Filtr metoda probkowania w dziedzinie czestotlwiosci
% filtr nierekursywny
rzad = 16;  % rzad filtra
b1 = fir1(rzad, wn1);
x_fir2 = filter(b1, 1, x);
% Obliczenie widma i mocy za pomoc¹ funkcji z lab 4.
[f_w, Moc, Wid] = fft_from_signal([x; x_fir2], fs);

% Wykresy
figure(3);
sgtitle(['Filtr dolnoprzepustowy FIR2, rzad=' num2str(rzad)...
    ' f_o=' num2str(fo1)]);
subplot(211);
plot(t, x); title('Sygna³ wymuszaj¹cy');
xlabel('Czas [s]'); ylabel('Amp.'); grid;
subplot(212);
plot(t, x_fir2); title('Sygna³ po filtracji');
xlabel('Czas [s]'); ylabel('Amp.'); grid;

figure(4);
sgtitle(['Filtr dolnoprzepustowy FIR2, rzad=' num2str(rzad)...
    ', f_o=' num2str(fo1)]);
subplot(211);
plot(f_w, Wid(1,:)); title('Widmo sygna³u wymuszaj¹cego');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;
subplot(212);
plot(f_w, Wid(2,:)); title('Widmo po filtracji');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;

% Wykresy porównawcze filtracji
figure(5);
spectrogram(x, 'yaxis'); title('Sygna³ wymuszaj¹cy');
figure(6);
spectrogram(x_filtered, 'yaxis'); title('Sygna³ po filtracji - IIR');
figure(7);
spectrogram(x_fir2, 'yaxis'); title('Sygna³ po filtracji - FIR');

% Filtrowanie i krótkoczasowa analiza czêstotliwoœciowa
% sygna³ów o zmiennej czestotliwoœci
T = 0:(1/fs):1.023;
% chirp(wektor czasu, czestestotliwosc start, koniec narastania,
% czestotliwoscdocelowa)
X = chirp(T,50,1.023,450);      % funkcj generuj¹ca sygna³
% filtracja przez filtry FIR i IIR
X_fir2 = filter(b1, 1, X);
X_iir = filter(B, A, X);

figure(8)
subplot(311); plot(T, X);
title('Sygna³ wymuszaj¹cy'); xlabel('Czas [s]'); ylabel('Amp.'); grid;

subplot(312); plot(T, X_fir2);
title('Filtracja FIR'); xlabel('Czas [s]'); ylabel('Amp.'); grid;

subplot(313); plot(T, X_iir);
title('Filtracja IIR'); xlabel('Czas [s]'); ylabel('Amp.'); grid;


figure(9)
% SFFT od sygnalu X, z oknem 256, zak³adk¹ okien 250, i 256 punktow¹ FFT
spectrogram(X, 256, 250, 256, 1E3, 'yaxis'); title('Sygna³ wymuszaj¹cy');

figure(10)
spectrogram(X_fir2, 256, 250, 256, 1E3, 'yaxis');
title('Sygna³ po fitlracji FIR');

figure(11)
spectrogram(X_iir, 256, 250, 256, 1E3, 'yaxis');
title('Sygna³ po fitlracji IIR');


%% Filtr œrodkowoprzepustowy
% Filtr Butterwortha
rzad1 = 8;          % rz¹d filtra
[B, A] = butter(rzad1, [wn2 wn1]);   % dla czestosci odciecia nr 1
x_filtered = filter(B, A, x);   % filtracja syg. przez uk³ad

% Filtr metoda probkowania w dziedzinie czestotlwiosci
% filtr nierekursywny
rzad2 = 16;  % rzad filtra
b1 = fir1(rzad2, [wn2 wn1]);
x_fir2 = filter(b1, 1, x);
% Obliczenie widma i mocy za pomoc¹ funkcji z lab 4.
[f_w, Moc, Wid] = fft_from_signal([x; x_filtered; x_fir2], fs);

% Wykresy
figure(12);
sgtitle(['Filtr srodkowoprzepustowy f_o=' num2str(fo1)...
     ', ' num2str(fo2) ']']);
subplot(311);
plot(t, x); title('Sygna³ wymuszaj¹cy');
xlabel('Czas [s]'); ylabel('Amp.'); grid;
subplot(312);
plot(t, x_filtered);
title(['Sygna³ po filtracji, Butterworth rzad=' num2str(rzad1)]);
xlabel('Czas [s]'); ylabel('Amp.'); grid;
subplot(313);
plot(t, x_fir2);
title(['Sygna³ po filtracji, FIR2 rzad=' num2str(rzad2)]);
xlabel('Czas [s]'); ylabel('Amp.'); grid;


figure(14);
sgtitle(['Filtr srodkowoprzepustowy f_o=' num2str(fo1)...
     ', ' num2str(fo2) ']']);
subplot(311);
plot(f_w, Wid(1,:)); title('Widmo sygna³u wymuszaj¹cego');
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;
subplot(312);
plot(f_w, Wid(2,:));
title(['Widmo po filtracji, Butterworth rzad=' num2str(rzad2)]);
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;
subplot(313);
plot(f_w, Wid(3,:));
title(['Widmo po filtracji, FIR2 rzad=' num2str(rzad2)]);
xlabel('Czestotliwosc [Hz]'); ylabel('Amp.'); grid;

% Wykresy porównawcze filtracji
figure(15);
spectrogram(x, 'yaxis'); title('Sygna³ wymuszaj¹cy');
figure(16);
spectrogram(x_filtered, 'yaxis'); title('Sygna³ po filtracji - IIR');
figure(17);
spectrogram(x_fir2, 'yaxis'); title('Sygna³ po filtracji - FIR');

% Filtrowanie i krótkoczasowa analiza czêstotliwoœciowa
% sygna³ów o zmiennej czestotliwoœci
T = 0:(1/fs):1.023;
% chirp(wektor czasu, czestestotliwosc start, koniec narastania,
% czestotliwoscdocelowa)
X = chirp(T,50,1.023,450);      % funkcj generuj¹ca sygna³
% filtracja przez filtry FIR i IIR
X_fir2 = filter(b1, 1, X);
X_iir = filter(B, A, X);

figure(18)
subplot(311); plot(T, X);
title('Sygna³ wymuszaj¹cy'); xlabel('Czas [s]'); ylabel('Amp.'); grid;

subplot(312); plot(T, X_fir2);
title('Filtracja FIR'); xlabel('Czas [s]'); ylabel('Amp.'); grid;

subplot(313); plot(T, X_iir);
title('Filtracja IIR'); xlabel('Czas [s]'); ylabel('Amp.'); grid;


figure(19)
% SFFT od sygnalu X, z oknem 256, zak³adk¹ okien 250, i 256 punktow¹ FFT
spectrogram(X, 256, 250, 256, 1E3, 'yaxis'); title('Sygna³ wymuszaj¹cy');

figure(20)
spectrogram(X_fir2, 256, 250, 256, 1E3, 'yaxis');
title('Sygna³ po fitlracji FIR');

figure(21)
spectrogram(X_iir, 256, 250, 256, 1E3, 'yaxis');
title('Sygna³ po fitlracji IIR');