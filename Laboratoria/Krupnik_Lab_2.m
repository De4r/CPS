%% Lab 2 - operacje wejscia wyjscia i parametry sygna³ów - Mateusz Krupnik
% Wyczyszczenie ekranu i generowanie przebiegów sinusoidalnych
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
% Zapis przebiegów do pliku, %12.4f - zapis wartosci o d³ 12 znaków, 4
% znaki precyzji, \n - nowy wiersz
uchwyt=fopen('dane1.txt','w');
fprintf(uchwyt,'%12.4f %12.4f %12.4f %12.4f %12.4f\n',[t;y1;y2;y3;y4]);
% Zamkniêcie pliku, i ponowne otwarcie, odczyt pliku do macierzy D
% %g - odczyt zapisu w postaci dziesietnej lub wykladniczje, usuniecie zer
% z konca zapisu, %e - notacja wykladnicza, %f - dzisiêtna
fclose('all');
uchwyt=fopen('dane1.txt','r');
D=fscanf(uchwyt,'%g %g %g %g %g \n',[5 inf]);
fclose('all');
% Zapis danych w postaci binarnej, a nastêpnie ich odczyranie, inf - odczyt
% do ostatniej kolumny
uchwyt1=fopen('dane1.bin','w');
fwrite(uchwyt1,[t;y1;y2;y3;y4],'float');
fclose('all');
uchwyt1=fopen('dane1.bin','r');
y5=fread(uchwyt1,[5 inf],'float');
fclose('all');

% Wykresy wygenerowanych przebiegów
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

%% Wynzaczanie parametrów sygna³ów za pomoc¹ stworzonych funkcji
% Wyznaczenie wartoœci minimalnej i maksymalnej
% Wywo³ywane funkcji signal_min i signal_max odpowiadaj¹ - min() i max()
a = [min(y1) min(y2) min(y3) min(y4)]
b = [max(y1) max(y2) max(y3) max(y4)]
% Wyznaczenie wartoœci œredniej za pomoc¹ funkcji signal_mean - mean()
y_mean = [mean(y1) mean(y3) mean(y3) mean(y4)]
% Energia sygnalu - za pomoc¹ funkcji signal_energy()
e = [signal_energy(y1) signal_energy(y2) ...
    signal_energy(y3) signal_energy(y4)]


function energy = signal_energy(signal, dt)
%	Signal energy
%	signal - signal, dt - time step
energy = 0;
if nargin > 1
    dt = dt;
else
    dt = 1;
end
for i=1:length(signal)
    energy = energy + (signal(i))^2*dt;
end

end

