% Lab 4. Tworzenie generatorow sygnalow podstawowych
% Wywolanie nastepuje blokowo.
% Do dzia�ania wymagane s� deklaracje funkcji:
% sbp, sup_1_2, swd, swj, tbp, tbpp, tup, tupp
% Deklaracje s� zakomentowane na ko�cu pliku gdyby 
% m pliki si� zgubi�y.
%                                   Mateusz Krupnik
clc; clear all; close all;
% Parametry
A=1;            % amplituda
f=1;            % czestotliwosc
fs=1000;        % czest. probkowania
t=0:(1/fs):10;  % wektor czasu
n=[7,30,99];    % wektor liczby probek
w=2*pi*f;       % wektor czestosci

%% Sygnal prostk�tny bipolarny
% Generowanie wykresu
y = sbp(w, A, t, n);

% Wykresy
figure(1)
sgtitle('Fala prostok�tna bipolarna');
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w) ...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal prost�ktny unipolarny wype�nienie 1/2
% Generowanie wykresu
y = sup_1_2(w, A, t, n);

% Wykresy
figure(2)
sgtitle('Fala prostok�tna unipolarna wyplenienie 1/2');
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w) ...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal prost�ktny unipolarny o dowolonym wypelnieniu
tau = 0.2; % okres, wypelnienie
% Generowanie wykresu
y = sup_wyp(f, A, t, n, tau);

% Wykresy
figure(3)
sgtitle(['Fala prostok�tna unipolarna wyplenienie ' num2str(tau)]);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w) ...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal trojkatny bipolarny
% Generowanie wykresu
y = tbp(w, A, t, n);

% Wykresy
figure(4)
sgtitle(['Fala trojkatna bipolarna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w)...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal trojkatny bipolarny piloksztaltny
% Generowanie wykresu
y = tbpp(w, A, t, n);

% Wykresy
figure(5)
sgtitle(['Fala trojkatna bipolarna pi�opkszta�tna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w) ...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal trojkatny unipolarny 
% Generowanie wykresu
y = tup(w, A, t, n);

% Wykresy
figure(6)
sgtitle(['Fala trojkatna unipolarna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w) ...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal trojkatny unipolarna piloksztaltna
% Generowanie wykresu
y = tupp(w, A, t, n);

% Wykresy
figure(7)
sgtitle(['Fala trojkatna unipolarna pi�okszta�tna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w) ...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal sinusoidalny wyprostowany dwupo��wkowy
% Generowanie wykresu
y = swd(w, A, t, n);

% Wykresy
figure(8)
sgtitle(['Fala sinusoidalna wyprostowana dwupo��wkowa']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w) ...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal sinusoidalny wyprostowany jednopo�owkowy
% Generowanie wykresu
y = swj(w, A, t, n);

% Wykresy
figure(9)
sgtitle(['Fala sinusoidalna wyprostowana jednopo�owkowa']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:));
    title(['Fala dla: w=' num2str(w) ...
        ' n=' num2str(n(i)) ' A=' num2str(A)]);
end


% %%%%%%%%%%% DEFINICJE FUNKCJI %%%%%%%%%%%%%%%%%%%

% Definicja funkcji 1
function x = sbp(w, A, t, n)
    % Funcja generuj�ca fale prostkatna bipolarna
    % w - czesto��, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t));
    for i=1:length(n)
        for j=1:2:2*n(i)
            x(i,:) = x(i,:) + ((1/j)*sin(j*w*t));
        end
    end
    x = x*4*A/pi;
end

% Definicja funkcji 2
function x = sup_1_2(w, A, t, n)
    % Funcja generuj�ca fale prostkatna unipolarna o wypelnieniu 1/2
    % w - cz�sto��, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t));
    for i=1:length(n)
        for j=1:4:2*n(i)
            x(i,:) = x(i,:) + ((1/j)*cos(j*w*t));
        end
        for j=3:4:2*n(i)
            x(i,:) = x(i,:) - ((1/j)*cos(j*w*t));
        end
    end
    x = x*2*A/pi + A/2;
end

% Definicja funkcji 3
function x = sup_wyp(f, A, t, n, tau)
    % Funcja generuj�ca fale prostk. unipolarna o dowolnym wypelnieniu
    % f - czestotliwosc, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t)); T = 1/f;
    for i=1:length(n)
        for j=1:n(i)
            x(i,:) = x(i,:) + ...
              sin(pi*j*tau/T)*cos(2*j*pi*f*t)/(pi*j*tau/T);
        end       
    end
    x = A*tau/T + 2*A*tau*x/T;
end

% Definicja funkcji 4
function x = tbp(w, A, t, n)
    % Funcja generuj�ca fale trojkatna bipolarna
    % w - cz�sto��, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t));
    for i=1:length(n)
        for j=1:4:n(i)
            x(i,:) = x(i,:) + ((1/j^2)*sin(j*w*t));
        end
        for j=3:4:n(i)
            x(i,:) = x(i,:) - ((1/j^2)*sin(j*w*t));
        end      
    end
    x = x*8*A/(pi^2);
end

% Definicja funkcji 5
function x = tbpp(w, A, t, n)
    % Funcja generuj�ca fale trojkatna bipolarna pilokszta�tna
    % w - cz�sto��, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t));
    for i=1:length(n)
        for j=1:2:n(i)
            x(i,:) = x(i,:) + ((1/j)*sin(j*w*t));
        end
        for j=2:2:n(i)
            x(i,:) = x(i,:) - ((1/j)*sin(j*w*t));
        end      
    end
    x = x*2*A/pi;
end

% Definicja funkcji 6
function x = tup(w, A, t, n)
    % Funcja generuj�ca fale trojkatna unipolarna
    % w - cz�sto��, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t));
    for i=1:length(n)
        for j=0:n(i)
            x(i,:) = x(i,:) + ((1/((2*j+1)^2))*cos((2*j+1)*w*t));
        end    
    end
    x = x*(-4*A)/(pi^2) + A/2;
end

% Definicja funkcji 7
function x = tupp(w, A, t, n)
    % Funcja generuj�ca fale trojkatna unipolarna pilokszta�tna
    % w - cz�sto��, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t));
    for i=1:length(n)
        for j=1:n(i)
            x(i,:) = x(i,:) - ((1/j)*sin(j*w*t));
        end   
    end
    x = x*A/pi + A/2;
end

% Definicja funkcji 8
function x = swd(w, A, t, n)
    % Funcja generuj�ca fale sinusoidalna wyprostowana dwupo�owkow�
    % w - cz�sto��, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t));
    for i=1:length(n)
        for j=1:n(i)
            x(i,:) = x(i,:) + (1/(4*j^2-1))*cos(2*j*w*t);
        end   
    end
    x = x*(-4)*A/pi + 2*A/pi;
end

% Definicja funkcji 9
function x = swj(w, A, t, n)
    % Funcja generuj�ca fale sinusoidalna wyprostowana jednopo�owkow�
    % w - cz�sto��, A - amplituda
    % t - wektor czasu, n - rzad ciagu
    x=zeros(length(n), length(t));
    for i=1:length(n)
        for j=1:n(i)
            x(i,:) = x(i,:) + (1/(4*j^2-1))*cos(2*j*w*t);
        end
        x(i,:) = x(i,:)*A*(-2)/pi + A/pi + sin(w*t)*A/2;
    end
end