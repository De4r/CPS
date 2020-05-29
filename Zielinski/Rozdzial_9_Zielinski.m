% Rodzial 9 Zielniski
%                       Mateusz Krupnik

% Implementacj algorytmu RADIX-2
clear all; close all; clc;
% Sygna³y
syg_test = 2        % 1 - 8 punktowy testowy, 2 - zdefinowany
% Sygna³ 1
if syg_test == 1
    N=8;                % liczba próbek sygna³u
    x=0:N-1;            % przyk³adowe wartoœci próbek
    f = 1:N;            % próbki
else
   N = 512;
   fpr = 1000; dt=1/fpr;
   t = dt*(0:N-1);
   x = sin(2*pi*150*t)+sin(2*pi*15*t);
   f = fpr*(0:N-1)/N;
end
xc = x;             % kopia sygna³u x

% obliczenie widma metod¹ zaimplemetowan¹ w Matlabie
wid_fft = fft(xc);


% Opcja 1 - odwracanie w implementacji 1 i FFT pogl¹dowe
syg = bitReverse1(x);       % decymacja próbek
wid_1 = fft_1(syg);         % Olbiczenie FFT
% Obliczenie b³êdów i wykresy
blad_real = abs(real(wid_1-wid_fft));
blad_imag = abs(imag(wid_1-wid_fft));
figure(1); subplot(311);
plot(f, abs(wid_1), f, abs(wid_fft), 'r'); title('Widmo wynikowe i orginalne');
legend('FFT Wlasne', 'FFT Matlab');
subplot(312); plot(blad_real); title('B³ad czêœci rzeczywsitej');
subplot(313); plot(blad_imag); title('B³ad czêœci urojonej');


%% Opcja 2 - odwracanie w implementacji 2 i FFT pogl¹dowe
syg = bitReverse2(x);       % decymacja próbek
wid_2 = fft_1(syg);         % Olbiczenie FFT
% Obliczenie b³êdów i wykresy
blad_real = abs(real(wid_2-wid_fft));
blad_imag = abs(imag(wid_2-wid_fft));
figure(2); subplot(311);
plot(f, abs(wid_1), f, abs(wid_fft), 'r'); title('Widmo wynikowe i orginalne');
legend('FFT Wlasne', 'FFT Matlab');
subplot(312); plot(blad_real); title('B³ad czêœci rzeczywsitej');
subplot(313); plot(blad_imag); title('B³ad czêœci urojonej');


%% Opcja 3 - odwracanie w implementacji 1 i FFT szybkie
syg = bitReverse1(x);       % decymacja próbek
wid_3 = fft_2(syg);         % Olbiczenie FFT
% Obliczenie b³êdów i wykresy
blad_real = abs(real(wid_3-wid_fft));
blad_imag = abs(imag(wid_3-wid_fft));
figure(3); subplot(311);
plot(f, abs(wid_1), f, abs(wid_fft), 'r'); title('Widmo wynikowe i orginalne');
legend('FFT Wlasne', 'FFT Matlab');
subplot(312); plot(blad_real); title('B³ad czêœci rzeczywsitej');
subplot(313); plot(blad_imag); title('B³ad czêœci urojonej');


%% Opcja 4 - odwracanie w implementacji 2 i FFT szybkie
syg = bitReverse2(x);       % decymacja próbek
wid_4 = fft_2(syg);         % Olbiczenie FFT
% Obliczenie b³êdów i wykresy
blad_real = abs(real(wid_4-wid_fft));
blad_imag = abs(imag(wid_4-wid_fft));
figure(4); subplot(311);
plot(f, abs(wid_1), f, abs(wid_fft), 'r'); title('Widmo wynikowe i orginalne');
legend('FFT Wlasne', 'FFT Matlab');
subplot(312); plot(blad_real); title('B³ad czêœci rzeczywsitej');
subplot(313); plot(blad_imag); title('B³ad czêœci urojonej');


%% Impletmentacje funkcji przestawiania próbek
function syg1=bitReverse1(syg_wej1)
    N1 = length(syg_wej1);  % liczba próbek
    MSB=log2(N1);	% liczba bitów numerów próbek
    for n=0:N1-1 	% kolejne próbki
        ncopy=n;	% stary numer próbki (kopia)
        nr=0;       % nowy numer próbki (inicjalizacja)
        for m=1:MSB % po wszystkich bitach
            if (rem(n,2)==0)	% czy jedynka na LSB
                n=n/2;          % jeœli nie, przesuñ w prawo
            else
                nr=nr+2^(MSB-m);	% dodaj 2^(MSB-m)
                n=(n-1)/2;      % odejmij jedynkê, przesuñ w prawo
            end
        end
        y(nr+1)=syg_wej1(ncopy+1);     % skopiuj we w³aœciwe miejsce
    end
    syg1 = y;          % podstaw wynik pod x
end

% Metoda przestawiana w miescu
function syg2=bitReverse2(syg_wej2)
    N2 = length(syg_wej2);  % liczba próbek
    a=1;
    for b=1:N2-1     % kolejne próbki
        if (b<a)	% porównanie czy indeks jest wiêkszy
            T=syg_wej2(a);
            syg_wej2(a)=syg_wej2(b);
            syg_wej2(b)=T;      % zamiana miejscami
        end
        c=N2/2;      % indeks po³owy próbek
        while (c<a)
            a=a-c; c=c/2;	% decymacja
        end
        a=a+c;
    end
    syg2=syg_wej2;
end

% Obliczanie FFT - wersja pogl¹dowa
% Wymaga podania sygna³u po decymacji funkcjami powy¿ej
function x=fft_1(x)
    Nsyg = length(x);    % liczba próbek
    for e = 1 : log2(Nsyg) % KOLEJNE ETAPY
        SM = 2^(e-1);	% szerokoœæ motylka
        LB = Nsyg/(2^e);	% liczba bloków
        LMB = 2^(e-1);	% liczba motylków w bloku
        OMB = 2^e;      % odleg³oœæ miêdzy blokami
        W = exp(-1i*2*pi/2^e);	% podstawa bazy Fouriera
        for b = 1 : LB          % KOLEJNE BLOKI
            for m = 1 : LMB     % KOLEJNE MOTYLKI
                g = (b-1)*OMB + m;	% indeks górnej próbki motylka
                d = g + SM;         % indeks dolnej próbki motylka
                xg = x(g);     % skopiowanie górnej próbki
                xd = x(d)*W^(m-1);	% korekta dolnej próbki
                % nowa górna próbka: górna plus dolna po korekcie
                x(g) = xg + xd;
                % nowa dolna próbka: górna minus dolna po korekcie
                x(d) = xg - xd; 
            end % koniec pêtli motylków
        end % koniec pêtli bloków
    end % koniec pêtli etapów
end

% Obliczanie FFT - wersja szybsza
% Wymaga podania sygna³u po decymacji funkcjami powy¿ej
function x=fft_2(x)
    Ns = length(x);    % liczba próbek
    for e = 1 : log2(Ns) % KOLEJNE ETAPY
        L = 2^e;        % d³ugoœæ bloków DFT, przesuniêcie bloków
        M = 2^(e-1);	% liczba motylków w bloku, szerokoœæ ka¿dego motylka
        Wi = 1;         % startowa wartoœæ wsp. bazy w etapie     
        W = cos(2*pi/L)-1i*sin(2*pi/L);	% mno¿nik bazy Fouriera
        for m = 1 : M           % KOLEJNE MOTYLKI
            for g = m : L : Ns	% W KOLEJNYCH BLOKACH
                d = g+M;        % g ? „górny”, d ? „dolny” indeks próbki motylka
                T2 = x(d)*Wi;	% „serce” FFT
                x(d) = x(g)-T2;	% nowa dolna próbka: górna minus „serce”
                x(g) = x(g)+T2;	% nowa górna próbka: górna plus „serce”
            end                 % koniec pêtli bloków
            Wi=Wi*W;            % kolejna wartoœæ bazy Fouriera
        end % koniec pêtli motylków
    end
end