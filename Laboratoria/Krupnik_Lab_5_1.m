% Lab 5. Procedura szybiej transfomraty Fouriera
% Defincija funkcji obliczaj¹cej FFT z sygna³ów
% znajduje sie na koncu pliku.
%                                   Mateusz Krupnik

% Generowanie sygna³u sinusoidalnego oraz jego widmo 
clc; clear all; close all;
f1=100;
fs=1000;
t=0:(1/fs):1.3;
Ts=1/fs;
A=1;
y=A*sin(2*pi*f1*t);
N=1024;
fft_moc=fft(y(1:N));
moc_wid=fft_moc.*conj(fft_moc)/N;
f=fs*(0:N/2-1)/N;
figure(1)
plot(f,moc_wid(1:N/2))

%% Analiza sygna³ów z Lab_4
% Dane podstawowe
clc; clear all; close all;
A=1; f=5; fs=1000; w=2*pi*f;
t=0:(1/fs):1;
n=[7,30,99];
j=1; % numer wykresu
%% Sygna³ prostok¹tny bipolarny
% Generowanie wykresu
y = sbp(w, A, t, n);

% Wykresy
figure(1);
sgtitle('Fala prostok¹tna bipolarna');
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czas [s]');
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(2);
sgtitle('Analiza widmowa fali prostok¹tnej bipolarnej');
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end

%% Sygna³ prostok¹tny unipolarny o wype³nieniu 1/2
% Generowanie wykresu
y = sup_1_2(w, A, t, n);

% Wykresy
figure(3);
sgtitle('Fala prostok¹tna unipolarna o wyp. 1/2');
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czas [s]');
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(4);
sgtitle('Analiza widmowa fali prostok¹tnej unipolarnej wyp. 1/2');
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end

%% Sygna³ prostok¹tny unipolarny o wype³nieniu dowolnym
% Generowanie wykresu
tau = 0.75; % okres, wypelnienie w procentach
if tau <=1 && tau >=0
    tau = tau*1/f; % w sekundach
else
    tau = 0.3
end
y = sup_wyp(f, A, t, n, tau);

% Wykresy
figure(5);
sgtitle(['Fala prostok¹tna unipolarna wyplenienie ' num2str(tau)]);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czas [s]');
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(6);
sgtitle(['Analiza widmowa fali prostok¹tnej unipolarnej o wyplenieniu ' num2str(tau)]);
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end

%% Sygnal trojkatny bipolarny
% Generowanie wykresu
y = tbp(w, A, t, n);

% Wykresy
figure(7)
sgtitle(['Fala trojkatna bipolarna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czas [s]');
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(8);
sgtitle(['Analiza widmowa fali trojkatnej bipolarnej']);
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end

%% Sygnal trojkatny bipolarny piloksztaltny
% Generowanie wykresu
y = tbpp(w, A, t, n);

% Wykresy
figure(9)
sgtitle(['Fala trojkatna bipolarna pi³opkszta³tna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(10);
sgtitle(['Analiza widmowa fali trojkatnej bipolarnej pi³opkszta³tnej']);
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end

%% Sygnal trojkatny unipolarny 
% Generowanie wykresu
y = tup(w, A, t, n);

% Wykresy
figure(11)
sgtitle(['Fala trojkatna unipolarna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(12);
sgtitle(['Analiza widmowa fali trojkatnej unipolarnej']);
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end

%% Sygnal trojkatny unipolarna piloksztaltna
% Generowanie wykresu
y = tupp(w, A, t, n);

% Wykresy
figure(13)
sgtitle(['Fala trojkatna unipolarna pi³okszta³tna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(14);
sgtitle(['Analiza widmowa fali trojkatnej pi³okszta³tnej']);
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end

%% Sygnal sinusoidalny wyprostowany dwupo³ówkowy
% Generowanie wykresu
y = swd(w, A, t, n);

% Wykresy
figure(15)
sgtitle(['Fala sinusoidalna wyprostowana dwupo³ówkowa']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(16);
sgtitle(['Analiza widmowa fali sinusoidalnej wyprostowanej dwupo³ówkowej']);
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end

%% Sygnal sinusoidalny wyprostowany jednopo³owkowy
% Generowanie wykresu
y = swj(w, A, t, n);

% Wykresy
figure(17)
sgtitle(['Fala sinusoidalna wyprostowana jednopo³owkowa']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

% Obliczenie FFT
[f_w, M, W] = fft_from_signal(y, fs);
figure(18);
sgtitle(['Analiza widmowa fali sinusoidalnej wyprostowanej jednopo³owkowej']);
for i=1:length(n)
    subplot(2,length(n),i)
    plot(f_w, M(i, :)); title(['Moc widmowa fali: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Moc widmowa');
    subplot(2,length(n),3+i)
    plot(f_w, W(i, :)); title(['Widmo dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
    xlabel('Czestotliwosc [Hz]'); ylabel('Widmo');
end


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