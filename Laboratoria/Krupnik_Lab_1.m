%% Lab 1 - operacje wejscia wyjscia - Mateusz Krupnik
% Generowanie przebiegu sinosuidalnego
clc; close all; clear all;
t=0:0.001:1;
A=0.7;
f=100;
omega=2*pi*f;
y=A*sin(omega*t);
% wykres 
figure(1)
plot(t,y)
% Otwarcie pliku w trybie zapisywania
uchwyt=fopen('uchwyt.txt','w');
% Zapis kolumny czasu i przebiegu sinusoidalnego 
fprintf(uchwyt,'%12.4f %12.4f\n',[t;y]);
% Zamkniêcie wsztstkich plików
fclose('all');
% Otworzenie pliku w trybie odczytu i wczytanie wartoœci do macierzy DANE
uchwyt=fopen('uchwyt.txt','r');
DANE=fscanf(uchwyt,'%g %g \n',[2 inf]);
fclose('all');
% Wykres
figure(2)
plot(DANE(1,:),DANE(2,:))