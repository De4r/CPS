%% Lab 3 - praca z plikiem dziwiêkowym - Mateusz Krupnik
clc; clear all; close all;
% Generowanie przebiegu, zapis i odczyt w postaci pliku wav
A=0.5;
B=-0.3;
f1=700;
fs=10000;
t=0:(1/fs):1;
y1=A*sin(2*pi*f1*t);
audiowrite('plik.wav', y1, fs);
clear all; % usuniecie danych
%% Odczyt danych
% odczyt danych z pliku wav i jego odtworzenie
[y, Fs] = audioread('plik.wav'); sound(y); pause(1);

% Odtworzenie osi czasu i generacja sygna³u
t = 0:(1/Fs):1;
f1 = 700; f2 = 70;
A = 1; y1 = A*sin(2*pi*f1*t); sound(y1); pause(t(end));
alfa1 = 2; alfa2 = 6;

%% T³umienie wyk³adnicze sygna³u
% Generacja sygna³u, odtworzenie dzwiêku i wykres
yt=y1.*exp(-alfa2*t);
sound(yt,Fs)
figure(1)
plot(t,yt);

%% Modulacja sygna³ów
% Generowanie sygna³ów zmodulowanych
ym=2*A*y1.*sin(2*pi*f2*t);
ym1=sin(2*pi*10*t).*sin(2*pi*1000*t);
ym2=sin(2*pi*10*t).*sin(2*pi*1000*t).*exp(-alfa1*t);
ym3=sin(2*pi*10*t).*sin(2*pi*1000*t).*exp(-alfa2*t);
% Generowanie wykresów
figure(2)
subplot(2,2,1); plot(t, ym); title('Modulacja amplitudy');
subplot(2,2,2); plot(t, ym1); title('Modulacja  amplitudy');
subplot(2,2,3); plot(t, ym2);
title('Modulacja amplitudy wraz z t³umieniem');
subplot(2,2,4); plot(t, ym3);
title('Modulacja amplitudy wraz z t³umieniem');
sgtitle('Sygna³y modulowane');
% Odtworzenie sygna³ów
sound(ym); pause(t(end)); sound(ym1); pause(t(end));
sound(ym2); pause(t(end)); sound(ym3);

%% Po³¹czenie sygna³ow w 2 kana³y i zapis pliku wav
% Po³¹czenie przebiegów i ich zapis
Y = [ym2; ym3]';
audiowrite('plik2.wav', Y, Fs);
% Odczyt z pliku, wykresy kana³ów oraz odtworzenie dzwiêku
[Y1, Fs] = audioread('plik2.wav');
figure(3)
plot(Y1(:, 1)); hold on; plot(Y1(:, 2), 'r');
hold off; title('Polaczone sygnaly');
sound(Y1, Fs); pause(t(end));

%% Po³¹czone sygna³y: sygna³ t³umiony i narastaj¹cy
% Generowanie przebiegu sygna³u narastaj¹cego oraz jego zapis
ym4 = y1.*(1-exp(-alfa1*t));
Y2 = [ym2; ym4]';
audiowrite('plik3.wav', Y2, Fs);
% Odczyt sygna³u, wykres i odtworzenie dzwiêku
[Y3, Fs] = audioread('plik3.wav');
figure(4)
plot(Y3(:, 1)); hold on; plot(Y3(:, 2), 'r.-');
hold off; title('Polaczone sygnaly');
sound(Y3, Fs); pause(t(end));