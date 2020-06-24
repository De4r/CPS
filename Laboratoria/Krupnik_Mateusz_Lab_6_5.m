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

%% Projekt filtra cyfrowego i jego dzialanie
% Parametry, niektóre z poprzednich filtrów
format long e;  % Typ formatowania zmiennych, 16 miejsc, wykladniczo
F_=[0 0.1 0.2 0.5 0.7 1];
M_=[1 1 1 0 0 0];
% Parametr rzêdu z poprzednich sekcji
[l, m] = yulewalk(N(1,rzad), F_, M_);     % Metoda Yule-Walkera
b = fir2(128, F_, M_); fir2(128, F_, M_); % Metoda próbkowania w d. czest.


[h, w] = freqz(b, 1);
Mag = 20*log10(abs(h)); Fi = phase(h)*180/pi; w = w/pi;
[h, w_] = freqz(l, m);
Mag_ = 20*log10(abs(h)); Fi_ = phase(h)*180/pi; w_ = w_/pi;

% Wykresy
figure(21);
sgtitle('Porownanie roznych metod projektowania');
subplot(211); plot(w, Mag, w_, Mag_);
legend('FIR2', 'Yule-Walker'); grid;
xlabel('Czest. znorm.'); ylabel('Amplituda [dB]');
title('Charakterystyka amplitudowa');
subplot(212); plot(w, Fi, w_, Fi_);
legend('FIR2', 'Yule-Walker'); grid;
xlabel('Czest. znorm.'); ylabel('Faza [\circ]');
title('Charakterystyka fazowa');