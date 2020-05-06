% Lab 4. Tworzenie generatorow sygnalow podstawowych
% Wywolanie nastepuje blokowo.
% Do dzia³ania wymagane s¹ deklaracje funkcji:
% sbp, sup_1_2, swd, swj, tbp, tbpp, tup, tupp
% Deklaracje s¹ zakomentowane na koñcu pliku gdyby 
% m pliki siê zgubi³y.
%                                   Mateusz Krupnik
clc; clear all; close all;
% Parametry
A=1;            % amplituda
f=1;            % czestotliwosc
fs=1000;        % czest. probkowania
t=0:(1/fs):10;  % wektor czasu
n=[7,30,99];    % wektor liczby probek
w=2*pi*f;       % wektor czestosci

%% Sygnal prostk¹tny bipolarny
% Generowanie wykresu
y = sbp(w, A, t, n);

% Wykresy
figure(1)
sgtitle('Fala prostok¹tna bipolarna');
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal prost¹ktny unipolarny wype³nienie 1/2
% Generowanie wykresu
y = sup_1_2(w, A, t, n);

% Wykresy
figure(2)
sgtitle('Fala prostok¹tna unipolarna wyplenienie 1/2');
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal prost¹ktny unipolarny o dowolonym wypelnieniu
tau = 0.2; % okres, wypelnienie
% Generowanie wykresu
y = sup_wyp(f, A, t, n, tau);

% Wykresy
figure(3)
sgtitle(['Fala prostok¹tna unipolarna wyplenienie ' num2str(tau)]);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal trojkatny bipolarny
% Generowanie wykresu
y = tbp(w, A, t, n);

% Wykresy
figure(4)
sgtitle(['Fala trojkatna bipolarna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal trojkatny bipolarny piloksztaltny
% Generowanie wykresu
y = tbpp(w, A, t, n);

% Wykresy
figure(5)
sgtitle(['Fala trojkatna bipolarna pi³opkszta³tna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal trojkatny unipolarny 
% Generowanie wykresu
y = tup(w, A, t, n);

% Wykresy
figure(6)
sgtitle(['Fala trojkatna unipolarna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal trojkatny unipolarna piloksztaltna
% Generowanie wykresu
y = tupp(w, A, t, n);

% Wykresy
figure(7)
sgtitle(['Fala trojkatna unipolarna pi³okszta³tna']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal sinusoidalny wyprostowany dwupo³ówkowy
% Generowanie wykresu
y = swd(w, A, t, n);

% Wykresy
figure(8)
sgtitle(['Fala sinusoidalna wyprostowana dwupo³ówkowa']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end

%% Sygnal sinusoidalny wyprostowany jednopo³owkowy
% Generowanie wykresu
y = swj(w, A, t, n);

% Wykresy
figure(7)
sgtitle(['Fala sinusoidalna wyprostowana jednopo³owkowa']);
for i=1:length(n) 
    subplot(length(n),1,i)
    plot(t, y(i,:)); title(['Fala dla: w=' num2str(w) ' n=' num2str(n(i)) ' A=' num2str(A)]);
end



% %%%%%%%%%%% DEFINICJE FUNKCJI %%%%%%%%%%%%%%%%%%%
% 
% 
% % Definicja funkcji 1
% function x = sbp(w, A, t, n)
%     % Funcja generuj¹ca fale prostkatna bipolarna
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t));
%     for i=1:length(n)
%         for j=1:2:n(i)
%             x(i,:) = x(i,:) + ((1/j)*sin(j*w*t));
%         end
%     end
%     x = x*4*A/pi;
% end
% 
% % Definicja funkcji 2
% function x = sup_1_2(w, A, t, n)
%     % Funcja generuj¹ca fale prostkatna unipolarna o wypelnieniu 1/2
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t));
%     for i=1:length(n)
%         for j=1:4:n(i)
%             x(i,:) = x(i,:) + ((1/j)*cos(j*w*t));
%         end
%         for j=3:4:n(i)
%             x(i,:) = x(i,:) - ((1/j)*cos(j*w*t));
%         end
%     end
%     x = x*2*A/pi + A/2;
% end
% 
% % Definicja funkcji 3
% function x = sup_wyp(f, A, t, n, tau)
%     % Funcja generuj¹ca fale prostkatna unipolarna o dowolnym wypelnieniu
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t)); T = 1/f;
%     for i=1:length(n)
%         for j=1:n(i)
%             x(i,:) = x(i,:) + sin(pi*j*tau/T)*cos(2*j*pi*f*t)/(pi*j*tau/T);
%         end       
%     end
%     x = A*tau/T + 2*A*tau*x/T;
% end
% 
% % Definicja funkcji 4
% function x = tbp(w, A, t, n)
%     % Funcja generuj¹ca fale trojkatna bipolarna
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t));
%     for i=1:length(n)
%         for j=1:4:n(i)
%             x(i,:) = x(i,:) + ((1/j^2)*sin(j*w*t));
%         end
%         for j=3:4:n(i)
%             x(i,:) = x(i,:) - ((1/j^2)*sin(j*w*t));
%         end      
%     end
%     x = x*8*A/(pi^2);
% end
% 
% % Definicja funkcji 5
% function x = tbpp(w, A, t, n)
%     % Funcja generuj¹ca fale trojkatna bipolarna pilokszta³tna
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t));
%     for i=1:length(n)
%         for j=1:2:n(i)
%             x(i,:) = x(i,:) + ((1/j)*sin(j*w*t));
%         end
%         for j=2:2:n(i)
%             x(i,:) = x(i,:) - ((1/j)*sin(j*w*t));
%         end      
%     end
%     x = x*2*A/pi;
% end
% 
% % Definicja funkcji 6
% function x = tup(w, A, t, n)
%     % Funcja generuj¹ca fale trojkatna unipolarna
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t));
%     for i=1:length(n)
%         for j=0:n(i)
%             x(i,:) = x(i,:) + ((1/((2*j+1)^2))*cos((2*j+1)*w*t));
%         end    
%     end
%     x = x*(-4*A)/(pi^2) + A/2;
% end
% 
% % Definicja funkcji 7
% function x = tupp(w, A, t, n)
%     % Funcja generuj¹ca fale trojkatna unipolarna pilokszta³tna
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t));
%     for i=1:length(n)
%         for j=1:n(i)
%             x(i,:) = x(i,:) + ((1/j)*sin(j*w*t));
%         end   
%     end
%     x = x*A/pi + A/2;
% end
% 
% % Definicja funkcji 8
% function x = swd(w, A, t, n)
%     % Funcja generuj¹ca fale sinusoidalna wyprostowana dwupo³owkow¹
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t));
%     for i=1:length(n)
%         for j=1:n(i)
%             x(i,:) = x(i,:) + (1/(4*j^2-1))*cos(2*j*w*t);
%         end   
%     end
%     x = x*(-4)*A/pi + 2*A/pi;
% end
% 
% % Definicja funkcji 9
% function x = swj(w, A, t, n)
%     % Funcja generuj¹ca fale sinusoidalna wyprostowana jednopo³owkow¹
%     % w - czestotliwosc, A - amplituda
%     % t - wektor czasu, n - rzad ciagu
%     x=zeros(length(n), length(t));
%     for i=1:length(n)
%         for j=1:n(i)
%             x(i,:) = x(i,:) + (1/(4*j^2-1))*cos(2*j*w*t);
%         end
%         x(i,:) = x(i,:)*A*(-2)/pi + A/pi + sin(w*t)*A/2;
%     end
% 
% end
