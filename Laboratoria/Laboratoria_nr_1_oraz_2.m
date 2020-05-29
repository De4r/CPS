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
% Zamkni�cie wsztstkich plik�w
fclose('all');
% Otworzenie pliku w trybie odczytu i wczytanie warto�ci do macierzy DANE
uchwyt=fopen('uchwyt.txt','r');
DANE=fscanf(uchwyt,'%g %g \n',[2 inf]);
fclose('all');
% Wykres
figure(2)
plot(DANE(1,:),DANE(2,:))

%% Lab 2 - operacje wejscia wyjscia i parametry sygna��w - Mateusz Krupnik
% Wyczyszczenie ekranu i generowanie przebieg�w sinusoidalnych
clc
clear all

A=0.5;
B=-0.3;
f1=700;
f2=1200;
fs=10000;
t=0:(1/fs):1;
y1=A*sin(2*pi*f1*t);
y2=B*sin(2*pi*f2*t);
y3=y1+y2;
y4=y1-y2;
sound(y1,fs); pause(t(end));
sound(y2,fs); pause(t(end));
sound(y3,fs); pause(t(end));
sound(y4,fs); pause(t(end));
% Zapis przebieg�w do pliku, %12.4f - zapis wartosci o d� 12 znak�w, 4
% znaki precyzji, \n - nowy wiersz
uchwyt=fopen('dane1.txt','w');
fprintf(uchwyt,'%12.4f %12.4f %12.4f %12.4f %12.4f\n',[t;y1;y2;y3;y4]);
% Zamkni�cie pliku, i ponowne otwarcie, odczyt pliku do macierzy D
% %g - odczyt zapisu w postaci dziesietnej lub wykladniczje, usuniecie zer
% z konca zapisu, %e - notacja wykladnicza, %f - dzisi�tna
fclose('all');
uchwyt=fopen('dane1.txt','r');
D=fscanf(uchwyt,'%g %g %g %g %g \n',[5 inf]);
fclose('all');
% Zapis danych w postaci binarnej, a nast�pnie ich odczyranie, inf - odczyt
% do ostatniej kolumny
uchwyt1=fopen('dane1.bin','w');
fwrite(uchwyt1,[t;y1;y2;y3;y4],'float');
fclose('all');
uchwyt1=fopen('dane1.bin','r');
y5=fread(uchwyt1,[5 inf],'float');
fclose('all');

% Wykresy wygenerowanych przebieg�w
figure(3)
subplot(2,2,1)
plot(t,y1); title('y1');
subplot(2,2,2)
plot(t,y2,'r'); title('y2');
subplot(2,2,3)
plot(t,y3,'k'); title('y3');
subplot(2,2,4)
plot(t,y4,'g'); title('y4');
sgtitle('Dane wygenerowane')

% Wykresy danych odczytanych z pliku .txt
figure(4)
subplot(2,2,1)
plot(D(1,:),D(2,:))
subplot(2,2,2)
plot(D(1,:),D(3,:),'r')
subplot(2,2,3)
plot(D(1,:),D(4,:),'k')
subplot(2,2,4)
plot(D(1,:),D(5,:),'g')
sgtitle('Dane odczytane z .txt')

% Wykresy danych odczytanych z pliku .bin - wykres y5
figure(5)
subplot(2,2,1)
plot(t,y1); title('y1');
subplot(2,2,2)
plot(t,y2,'r'); title('y2');
subplot(2,2,3)
plot(t,y3,'k'); title('y3');
subplot(2,2,4)
plot(t,y4,'g'); title('Odczytana kolumna y4 z pliku .bin');
sgtitle('Dane odczytane z .bin')

%% Wynzaczanie parametr�w sygna��w za pomoc� stworzonych funkcji
% Wyznaczenie warto�ci minimalnej i maksymalnej
% Wywo�ywane funkcji signal_min i signal_max odpowiadaj� - min() i max()
a = [signal_min(y1) signal_min(y2) signal_min(y3) signal_min(y4)]
b = [signal_max(y1) signal_max(y2) signal_max(y3) signal_max(y4)]
% Wyznaczenie warto�ci �redniej za pomoc� funkcji signal_mean - mean()
y_mean = [signal_mean(y1) signal_mean(y3) ...
    signal_mean(y3) signal_mean(y4)]
% Energia sygnalu - za pomoc� funkcji signal_energy()
e = [signal_energy(y1) signal_energy(y2) ...
    signal_energy(y3) signal_energy(y4)]




