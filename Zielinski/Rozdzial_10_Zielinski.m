% Rodzial 10 Zielniski
%                       Mateusz Krupnik

%% Filtracja cyfrowa z wykorzystaniem bufor�w
%y = filter(b, a ,x);    % Funkcja wbudowana Matlaba
% filterBK i filterBP - definicje ponizej

%% Projektowanie rekursywnych filtr�w cyfrowych metod� zer i biegun�w
clear all; clc; close all;
% Parametry do rysowania okregu
NP = 1000; fi = 2 * pi * (0:1:NP-1)/NP; s = sin(fi); c = cos(fi);
fpr = 1000;     % czestotliwo�� probkowania
% syngaly testowe
Nx = 1024; n=0:Nx-1; dt=1/fpr; t=dt*n;
f1 = 10; f2 = 50; f3 = 250;
x = sin(2*pi*f1*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t);
xd = zeros(1, Nx); xd(1) = 1;
for i=1:2
   if (i == 1)      % FILTR LP
       filtr = "dolnoprzepustowy";
       fz = [ 50 ]; % Czestotliwo�� zer
       fp = [ 10 ]; % Czestotliwosc biegunow
       Rz = [ 1 ];	% wspoczynniki (promienie) zer
       Rp = [ 0.98 ];           % wspoczynniki (promienie) biegn�w (!= 1)
       fmax = 100; df = 0.1;    % widmo fouriera
   else
       filtr = "srodkoprzepustowy";
       fz = [ 50 100 150 350 400 450 ];	% Czestotliwo�� zer
       fp = [ 200 250 300 ];	% Czestotliwosc biegunow
       Rz = [ 1 ];              % wspoczynniki (promienie) zer
       Rp = [ 0.96 ];           % wspoczynniki (promienie) biegn�w (!= 1)
       fmax = 500; df = 0.1;    % widmo fouriera
   end
   
   % Obliczenie zer i biegun�w trnasmitancji H(z)
   fi_z = 2*pi*(fz/fpr);    % k�ty zer
   fi_p = 2*pi*(fp/fpr);    % k�ty biegunow
   z = Rz .* exp(1i*fi_z);  % zera
   p = Rp .* exp(1i*fi_p);  % bieguny
   z = [ z conj(z)]; p = [ p conj(p)];  % dodanie par sprz�onych
   
   % Po�o�eine zer i beign�w
   figure(6-*i - 5);
   plot(s, c, '-k', real(z), imag(z), 'or', real(p), imag(p), 'xb');
   title(["Zera i biguny filtra: " + filtr ]); legend('Zera','Bieguny');
   grid on;
   
   % Obliczenia wspo�czynnik�w a i b z zer i biegnow
   wzm = 1; [b, a] = zp2tf(z', p', wzm);
   
   % Charakterystyka czestotliwo�ciowa H(f)
   % czestotliwo��, czesto��, cz�sto�� unormowana
   f = 0 : df : fmax; w = 2*pi*f; wn = 2*pi*f/fpr;
   H = freqz(b, a, wn);
   Habs = abs(H); Hdb = 20*log10(Habs); Hfa = unwrap(angle(H));
   figure(6*i-4);
   sgtitle(["Analiza czestot. dla filtra: " + filtr]);
   subplot(311); plot(f, Habs); grid on; title('|H(f)|');
   xlabel('Czestotliwo�� f [Hz]'); ylabel('Ampl.');
   subplot(312); plot(f, Hdb); grid on; title('|H(f)| [dB]');
   xlabel('Czestotliwo�� f [Hz]'); ylabel('Ampl. [dB]');
   subplot(312); plot(f, Hfa); grid on; title('angle(H(f))');
   xlabel('Czestotliwo�� f [Hz]'); ylabel('rad.');
   
   % Filtracja synga�u x
   y1 = filter(b, a, x);
   y2 = filterBP(b, a, x);
   y3 = filterBK(b, a, x);
   
   figure(6*i-3)
   sgtitle(["Filtracja sygna�u sinusoidalnego dla filtra: " + filtr]);
   subplot(211); plot(t, x);
   title('Sygna� wej�ciowy x'); xlabel('Czas [s]');
   ylabel('Ampl.'); grid on;
   subplot(212); plot(t, y1, t, y2, t, y3);
   title('Sygna� wyj�ciowy y'); xlabel('Czas [s]');
   ylabel('Ampl.'); grid on;
   legend('Matlab filter', 'filterBP', 'filterBK');
   
   % Filtracja synga�u xd
   y1d = filter(b, a, xd);
   y2d = filterBP(b, a, xd);
   y3d = filterBK(b, a, xd);
   
   figure(6*i-2)
   sgtitle(["Filtracja sygna�u jednostkowego dla filtra: " + filtr]);
   subplot(211); plot(t, xd);
   title('Sygna� wej�ciowy xd'); xlabel('Czas [s]');
   ylabel('Ampl.'); grid on;
   subplot(212); plot(t, y1d, t, y2d, t, y3d);
   title('Sygna� wyj�ciowy y'); xlabel('Czas [s]');
   ylabel('Ampl.'); grid on;
   legend('Matlab filter', 'filterBP', 'filterBK');
   
   
   % Sygnaly w dziedzinie czestotliwosci
    n=Nx/2+1:Nx; X = freqz(x(n),1,wn)/(Nx/4);
    Y = freqz(y1(n),1,wn)/(Nx/4);
    Y(2, :) = freqz(y2(n),1,wn)/(Nx/4);
    Y(3, :) = freqz(y3(n),1,wn)/(Nx/4);
    X = abs(X); Y = abs(Y);
    
    figure(6*i-1); (["Filtracja sygna�u sinusoidalnego dla filtra: " + filtr]);
    subplot(211); plot(f,X);
    title('Wejscie X(f)'); ylabel('Ampl.'); grid on;
    subplot(212); plot(f, Y(1, :), f, Y(2, :), f, Y(3, :));
    title('Wyj�cie Y(f)'); ylabel('Ampl.'); grid on;
    xlabel('f [Hz]'); legend('Matlab filter', 'filterBP', 'filterBK');
    
    % Sygnal impulsowy
    
    Xd = freqz(xd(n),1,wn)/(Nx/4);
    Yd = freqz(y1d(n),1,wn)/(Nx/4);
    Yd(2, :) = freqz(y2d(n),1,wn)/(Nx/4);
    Yd(3, :) = freqz(y3d(n),1,wn)/(Nx/4);
    Xd = abs(Xd); Yd = abs(Yd);
    figure(6*i); (["Filtracja sygna�u impulsowego dla filtra: " + filtr]);
    subplot(211); plot(f, Xd);
    title('Wejscie X(f)'); ylabel('Ampl.'); grid on;
    subplot(212); plot(f, Yd(1, :), f, Yd(2, :), f, Yd(3, :));
    title('Wyj�cie Y(f)'); ylabel('Ampl.'); grid on;
    xlabel('f [Hz]'); legend('Matlab filter', 'filterBP', 'filterBK');
end









%%%%%%% DEFINICJE FUNKCJI %%%%%%%%%%%%%%%%%%

% Filtracj z buforami przesuwnymi
function y = filterBP(b, a, x)
	Nx = length(x);                 % D�ugo�� sygna�u x
    M = length(b); N = length(a);	% Dlugo�� bufor�w
    a = a(2:N); N = N - 1;          % Usuni�cie a_1 = 1
    bx = zeros(1, M); by = zeros(1, N); % Prealokacja bufor�w
    y = [];
    for n=1:Nx
        bx = [x(n) bx(1:M-1)];
        y(n) = sum(bx .* b) - sum(by .* a);
        by = [y(n) by(1:N-1)];
    end
end

% Filtracja z buforami kolowymi
function y = filterBK(b, a, x)
    Nx = length(x);                 % D�ugo�� sygna�u x
    M = length(b); N = length(a);	% Dlugo�� bufor�w
    a = a(2:N); N = N - 1;          % Usuni�cie a_1 = 1
    bx = zeros(1, M); by = zeros(1, N); % Prealokacja bufor�w
    y = [];
    ix = 1; iy = 1;                 % wska�niki od 1 miejsca tablic
    for n = 1 : Nx
        bx(ix) = x(n);              % Bufor wej�cia
        sum = 0; ib = 1; ia = 1;    % wskazniki
        
        for k = 1 : M - 1           % Sumowanie probek wejscia
            sum = sum +bx(ix)*b(ib);
            ix = ix - 1;
            if (ix == 0) ix = M; end
            ib = ib + 1;
        end
        sum = sum + bx(ix)*b(ib);   % ostatni skladnik sumy
        
        for k = 1 : N - 1
            sum = sum - by(iy)*a(ia);
            iy = iy - 1;
            if (iy == 0) iy = N; end
            ia = ia + 1;
        end
        sum = sum - by(iy)*a(ia);
        
        y(n) = sum;                 % wynik filtracji
        by(iy) = sum;               % dodanie wyj�cia do bufora wyj��
    end
end