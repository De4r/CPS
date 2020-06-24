% Rodzial 9 Zielinski
%                       Mateusz Krupnik

% Implementacj algorytmu RADIX-2
clear all; close all; clc;
% Sygna�y
syg_test = 2        % 1 - 8 punktowy testowy, 2 - zdefinowany
% Sygna� 1
if syg_test == 1
    N=8;                % liczba pr�bek sygna�u
    x=0:N-1;            % przyk�adowe warto�ci pr�bek
    f = 1:N;            % pr�bki
else
   N = 512;
   fpr = 1000; dt=1/fpr;
   t = dt*(0:N-1);
   x = sin(2*pi*150*t)+sin(2*pi*15*t);
   f = fpr*(0:N-1)/N;
end
xc = x;             % kopia sygna�u x

% obliczenie widma metod� zaimplemetowan� w Matlabie
wid_fft = fft(xc);


% Opcja 1 - odwracanie w implementacji 1 i FFT pogl�dowe
syg = bitReverse1(x);       % decymacja pr�bek
wid_1 = fft_1(syg);         % Olbiczenie FFT
% Obliczenie b��d�w i wykresy
blad_real = abs(real(wid_1-wid_fft));
blad_imag = abs(imag(wid_1-wid_fft));
figure(1); subplot(311);
plot(f, abs(wid_1), f, abs(wid_fft), 'r');
title('Widmo wynikowe i orginalne');
legend('FFT Wlasne', 'FFT Matlab');
subplot(312); plot(blad_real); title('B�ad cz�ci rzeczywsitej');
subplot(313); plot(blad_imag); title('B�ad cz�ci urojonej');


%% Opcja 2 - odwracanie w implementacji 2 i FFT pogl�dowe
syg = bitReverse2(x);       % decymacja pr�bek
wid_2 = fft_1(syg);         % Olbiczenie FFT
% Obliczenie b��d�w i wykresy
blad_real = abs(real(wid_2-wid_fft));
blad_imag = abs(imag(wid_2-wid_fft));
figure(2); subplot(311);
plot(f, abs(wid_1), f, abs(wid_fft), 'r');
title('Widmo wynikowe i orginalne');
legend('FFT Wlasne', 'FFT Matlab');
subplot(312); plot(blad_real); title('B�ad cz�ci rzeczywsitej');
subplot(313); plot(blad_imag); title('B�ad cz�ci urojonej');


%% Opcja 3 - odwracanie w implementacji 1 i FFT szybkie
syg = bitReverse1(x);       % decymacja pr�bek
wid_3 = fft_2(syg);         % Olbiczenie FFT
% Obliczenie b��d�w i wykresy
blad_real = abs(real(wid_3-wid_fft));
blad_imag = abs(imag(wid_3-wid_fft));
figure(3); subplot(311);
plot(f, abs(wid_1), f, abs(wid_fft), 'r');
title('Widmo wynikowe i orginalne');
legend('FFT Wlasne', 'FFT Matlab');
subplot(312); plot(blad_real); title('B�ad cz�ci rzeczywsitej');
subplot(313); plot(blad_imag); title('B�ad cz�ci urojonej');


%% Opcja 4 - odwracanie w implementacji 2 i FFT szybkie
syg = bitReverse2(x);       % decymacja pr�bek
wid_4 = fft_2(syg);         % Olbiczenie FFT
% Obliczenie b��d�w i wykresy
blad_real = abs(real(wid_4-wid_fft));
blad_imag = abs(imag(wid_4-wid_fft));
figure(4); subplot(311);
plot(f, abs(wid_1), f, abs(wid_fft), 'r');
title('Widmo wynikowe i orginalne');
legend('FFT Wlasne', 'FFT Matlab');
subplot(312); plot(blad_real); title('B�ad cz�ci rzeczywsitej');
subplot(313); plot(blad_imag); title('B�ad cz�ci urojonej');


%% Impletmentacje funkcji przestawiania pr�bek
function syg1=bitReverse1(syg_wej1)
    N1 = length(syg_wej1);  % liczba pr�bek
    MSB=log2(N1);	% liczba bit�w numer�w pr�bek
    for n=0:N1-1 	% kolejne pr�bki
        ncopy=n;	% stary numer pr�bki (kopia)
        nr=0;       % nowy numer pr�bki (inicjalizacja)
        for m=1:MSB % po wszystkich bitach
            if (rem(n,2)==0)	% czy jedynka na LSB
                n=n/2;          % je�li nie, przesu� w prawo
            else
                nr=nr+2^(MSB-m);	% dodaj 2^(MSB-m)
                n=(n-1)/2;      % odejmij jedynk�, przesu� w prawo
            end
        end
        y(nr+1)=syg_wej1(ncopy+1);     % skopiuj we w�a�ciwe miejsce
    end
    syg1 = y;          % podstaw wynik pod x
end

% Metoda przestawiana w miescu
function syg2=bitReverse2(syg_wej2)
    N2 = length(syg_wej2);  % liczba pr�bek
    a=1;
    for b=1:N2-1     % kolejne pr�bki
        if (b<a)	% por�wnanie czy indeks jest wi�kszy
            T=syg_wej2(a);
            syg_wej2(a)=syg_wej2(b);
            syg_wej2(b)=T;      % zamiana miejscami
        end
        c=N2/2;      % indeks po�owy pr�bek
        while (c<a)
            a=a-c; c=c/2;	% decymacja
        end
        a=a+c;
    end
    syg2=syg_wej2;
end

% Obliczanie FFT - wersja pogl�dowa
% Wymaga podania sygna�u po decymacji funkcjami powy�ej
function x=fft_1(x)
    Nsyg = length(x);    % liczba pr�bek
    for e = 1 : log2(Nsyg) % KOLEJNE ETAPY
        SM = 2^(e-1);	% szeroko�� motylka
        LB = Nsyg/(2^e);	% liczba blok�w
        LMB = 2^(e-1);	% liczba motylk�w w bloku
        OMB = 2^e;      % odleg�o�� mi�dzy blokami
        W = exp(-1i*2*pi/2^e);	% podstawa bazy Fouriera
        for b = 1 : LB          % KOLEJNE BLOKI
            for m = 1 : LMB     % KOLEJNE MOTYLKI
                g = (b-1)*OMB + m;	% indeks g�rnej pr�bki motylka
                d = g + SM;         % indeks dolnej pr�bki motylka
                xg = x(g);     % skopiowanie g�rnej pr�bki
                xd = x(d)*W^(m-1);	% korekta dolnej pr�bki
                % nowa g�rna pr�bka: g�rna plus dolna po korekcie
                x(g) = xg + xd;
                % nowa dolna pr�bka: g�rna minus dolna po korekcie
                x(d) = xg - xd; 
            end % koniec p�tli motylk�w
        end % koniec p�tli blok�w
    end % koniec p�tli etap�w
end

% Obliczanie FFT - wersja szybsza
% Wymaga podania sygna�u po decymacji funkcjami powy�ej
function x=fft_2(x)
    Ns = length(x);    % liczba pr�bek
    for e = 1 : log2(Ns) % KOLEJNE ETAPY
        L = 2^e;        % d�ugo�� blok�w DFT, przesuni�cie blok�w
        M = 2^(e-1);	% lba motylk�w w bloku, szeroko�� ka�dego motylka
        Wi = 1;         % startowa warto�� wsp. bazy w etapie     
        W = cos(2*pi/L)-1i*sin(2*pi/L);	% mno�nik bazy Fouriera
        for m = 1 : M           % KOLEJNE MOTYLKI
            for g = m : L : Ns	% W KOLEJNYCH BLOKACH
                % g �g�rny�, d �dolny� indeks pr�bki motylka
                d = g+M;        
                T2 = x(d)*Wi;	% �serce� FFT
                x(d) = x(g)-T2;	% nowa dolna pr�bka: g�rna minus �serce�
                x(g) = x(g)+T2;	% nowa g�rna pr�bka: g�rna plus �serce�
            end                 % koniec p�tli blok�w
            Wi=Wi*W;            % kolejna warto�� bazy Fouriera
        end % koniec p�tli motylk�w
    end
end