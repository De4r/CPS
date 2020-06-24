% Lab 6. 
% Dla filtra FIR sprawdzic jego dzia³anie.
%                                   Mateusz Krupnik
clear all; close all; clc;
%% Sprawdzenie dzia³ania filtrów cyfrowych
% Dane sygna³ów
f1=100; f2=250; f3=400; fs=1000;    % Czestotliwosci skladowych i Nyquista
A1=1; A2=0.8; A3=0.65;              % Aplitudy skladowych
t=0:(1/fs):1.023;                   % Wektor czasu
% Sygna³ wymuszenia - sk³adowa 3 harmonicznych
x=A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t)+A3*sin(2*pi*f3*t);

% Parametry
fo1 = 300; fo2 = 200;           % Czestotliwosci ociecia
wn1 = fo1*2/fs; wn2 = fo2*2/fs; % Znormalizowane czestotliwosci
wn = [wn1, wn2]; fo = [fo1, fo2];
% Filtr Butter - dolnoprzepustowy
N = [2, 4 ,8]; rzad = 3;    % wybór rzêdu filtra

%% Filtr srodkowoprzepustowy Butterwortha 
% Parametry z poprzedniego filtru
[l, m] = butter(N(1,rzad), fliplr(wn)); % Odwrócenie kolejnoœci cz. odc.
y = filter(l, m, x);                    % Odpowiedz filtra
[f_w, Moc, Wid] = fft_from_signal([x; y], fs);   % Funkcja z Lab 4 i 5

% Dzwiek
sound(x); pause(t(end)); sound(y); pause(t(end));
% Wykresy
figure(19);
sgtitle(['Filtr srodkowoprzepustowy f_{o1}=' ...
    num2str(fo(2)) ' f_{o2}=' num2str(fo(1))]);
subplot(321); plot(t, x);
xlabel('Czas [s]'); ylabel('Amplituda'); grid;
title('Sygna³ wymuszaj¹cy');
subplot(322); plot(t, y);
xlabel('Czas [s]'); ylabel('Amplituda'); grid;
title('Sygna³ po filtracji');
subplot(323); plot(f_w, Wid(1,:));
xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Widmo wymuszenia');
subplot(324); plot(f_w, Wid(2,:));
xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Widmo odpowiedzi');
subplot(325); plot(f_w, Moc(1,:));
xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Moc widmowa wymuszenia');
subplot(326); plot(f_w, Moc(2,:));
xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Moc widmowa odpowiedzi');

% DEFINICJA FUNKCI %%%%%%%%%%%%%%%%%
function [f, M, W] = fft_from_signal(y, fs)
%	fft_from_signal 
%   Summary of this function goes here
%   Detailed explanation goes here
%   y - signal matrix, with signals as rows
%   f - frequency, M - power sepctrum, W - spectrum
N = length(y);
for i=1:size(y, 1)
    fft_moc=fft(y(i, 1:N));
    moc_wid=fft_moc.*conj(fft_moc)/N;
    widmo=sqrt(fft_moc.*conj(fft_moc))/N;
    f=fs*(0:N/2-1)/N;
    M(i,:)=moc_wid;
    W(i,:)=widmo;
end
M = M(:, floor(1:N/2));
W = W(:, floor(1:N/2));
end