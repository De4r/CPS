% Rodzial 12 Zielniski
%                       Mateusz Krupnik

clc; close all; clear all;
% Zadanie 1
% Æwiczenie: Projektowanie nierekursywnych filtrów
% cyfrowych metod¹ próbkowania w dziedzinie czêstotliwoœci

N = 41; % dlugosc filtra
for i=1:4
    if (i==1 | i==3)
        M = (N -1) / 2;
        M2 = M/2;
        M4 = M/4;
        nn = N;
        Z4 = zeros(1, M4); J4 = ones(1, M4);
    else
        M = (N -1) / 2;
        M2 = M / 2;
        nn = 2 * M;
    end

    Z2 = zeros(1, M2); J2 = ones(1, M2);

    if i==1
        Ar = [ 1 J2 Z2 Z2 J2];
        Ai = zeros(1, nn);
        typ = "I";
    elseif i==2
        Ar = [J2 Z2 Z2 J2];
        Ai = zeros(1, nn);
        typ = "II";
    elseif i==3
        Ar = zeros(1, nn);
        Ai = [0 Z4 J2 Z4 Z4 -J2 Z4];
        typ = "III";
    else
        Ar = zeros(1, nn);
        Ai = [Z4 J2 Z4 Z4 -J2 Z4];
        typ = "IV";
    end
    A = Ar + 1i*Ai;
    
    n = 0 : nn-1; f = n/nn; h = zeros(1, nn);
    for k=0:nn-1
        h = h + A(k+1)*exp(1i*2*pi*k/nn*(n-M));
    end
    h = real(h/nn);
    ho = h.*blackman(nn)';
    
    % Wykresy filtrów
    figure(2*i-1)
    title(["Odp. impulsowa " + typ]);
    stem(n, h); hold on; stem(n, ho); hold off;
    legend('Okno prost.', 'Okno Blackmana');
    
    NF = 500; k=0:NF-1; fn=k/NF; wn=2*pi*fn;
    for k=0:NF-1
        temp = exp(-1i*2*pi*k/NF*(n-M))
        H(k+1) = temp * h';
        Ho(k+1) = temp * ho';
    end
    figure(i);
    if (i==1 | i==2) Ax = Ar; Hx=real(H); Hxo=real(Ho); end % dla filtra typu I i II
    if (i==3 | i==4) Ax = Ai; Hx=imag(H); Hxo=imag(Ho); end % dla filtra typu III i IV
    figure(2*i)
    sgtitle("Charakterystyki filtra typu " + typ);
    subplot(211);
    plot(f, Ax, 'ob', fn, Hx, fn, Hxo);
    grid; title('real(H) lub imag(H)');
    legend('Wymagania','Okno prost.', 'Okno Blackmana');
    subplot(212);
    plot(fn,20*log10(abs(H)), fn, 20*log10(abs(Ho)));
    grid; title('Modu³ |H| w dB');
    xlabel('Czestot. [Hz]'); ylabel('Ampl.');
    legend('Okno prost.', 'Okno Blackmana');
end

%% Æwiczenie: Projektowanie nierekursywnych filtrów cyfrowych
% metod¹ próbkowania w dziedzinie
% czêstotliwoœci i optymalizacji œredniokwadratowej
% Projektowanie filtrów FIR metod¹ WA¯ONEJ minimalizacji b³êdu œredniokwadratowego
% pomiêdzy zadan¹ charakterystyk¹, spróbkowan¹ w dziedzinie czêstotliwoœci, a charakterystyk¹ otrzymywan¹

% Wymagania dla filtra cyfrowego
M = 20; % po³owa d³ugoœci filtra N = 2*M + 1
K = 50; % liczba punktow charakterystyki >2M
% Wymagania czestotliwoœciowe, Ak
L1 = floor(K/4);    % podzia³ na 4 
% jedynki, przesjsciowe punkty, zera, przejsciowe punkty, jedynki
Ak = [ ones(1, L1) 0.6 0.4 zeros(1, K-(2*L1-1)-4) 0.4 0.6 ones(1, L1-1)]';

% Wagi punktow charakterystyki
wp = 1; % waga odcinka Pass
wt = 1; % waga odcinka przejsciowego
ws = 1; % waga odcinka Stop
% Wektor wag
w = [wp*ones(1, L1) wt wt ws*ones(1, K-(2*L1-1)-4) wt wt wp*ones(1, L1-1)];
W = diag(w);

% Wyznacznie macierzy F -> W*F*h=W*(Ak + e)
F = [];
n = 0:M-1;
for k = 0 : K - 1
    F = [F; 2*cos(2*pi*(M-n)*k/K) 1];
end
% wyznaczenie odp impulsowej
h = (W*F) \ (W*Ak);
h = [ h; h(M:-1:1) ];
% Wykresy
n = 0:2*M;
figure(57);
stem(n, h); grid on;
title('Odpowiedz impulsowa filtra'); xlabel('Probki n'); ylabel('Ampl.');
NF = 500; wn = 0:pi/(NF-1):pi; fn = wn/(2*pi);
H = freqz(h, 1, wn);
figure(58);
sgtitle('Charakterystyka filtra cyfrowego metod¹ wa¿onej minimalizacji b³edu');
subplot(311); plot(fn, abs(H)); grid; title('Modu³ odpowiedzi czest.');
xlabel('Czestotliwosc znorm. [Hz]'); ylabel('Ampl.');
subplot(312); plot(fn, 180/pi*unwrap(angle(H))); grid; title('Faza odpowiedzi czest.');
xlabel('Czestotliwosc znorm. [Hz]'); ylabel('deg.');
subplot(313); plot(fn, 20*log10(abs(H))); grid; title('Modul odpowiedzi czest.');
xlabel('Czestotliwosc znorm. [Hz]'); ylabel('Ampl. [dB]');

%% Æwiczenie: Projektowanie nierekursywnych filtrów cyfrowych w
% dziedzinie czêstotliwoœci metod¹
% aproksymacji Czebyszewa (algorytm Remeza)

L = 20; % liczba wspolczynnikow N=2L-1
Nr = 5; % szerokoœæ pasma przepustowego 0 < Nr < L
wp = 1; ws = 1; %wagi pasm Pass i Stop
R = 200;    % Zbior testowy/zbio ekstremow
tol = 10^-8;    % toleracnja bledu

% Parametry
M = L + 1;  % liczba czestotliwosci ekstremow
K = 1+R*(M-1);  % Liczba badanych czestotliwosci
fz = (0 : K-1)/(K-1);   % zbior czestotliwosci znorm.
k11 = 1; k12 = 1+Nr*R;  % granice pasma Pass
k21 = 1+(Nr+1)*R; k22 = K;  % granice pasma Stop
K1 = 1+Nr*R+R/2;            % nr probki dla czestot. granicznej
fd = [fz(1:K1) fz(K1:K)];   % czest. charakter filtra
Hd = [ones(1, K1) zeros(1, K-K1)];  % wzmocnienia
Wg = [wp*ones(1, K1) ws*ones(1, K-K1)]; % wagi
i_maximum = 1:R:K;      % inkdesy czestotliwosci ekstremow

% Wybor startowego zbioru czestotliwoœci ekstremów
feMAX = fz(1:R:K);
sigmaMAX = 10^15; sigma = 0;

% obliczenia
n = 0 : L-1;
while(sigmaMAX-sigma > tol)
    disp(['Zmniejszenie bledu: ' num2str(sigmaMAX-sigma)]);
    H = Hd(i_maximum)'; W = Wg(i_maximum);
    fe = feMAX;
    % macierz cosinusów
    A = [];
    for m = 0:M-1
        A = [A; cos(pi*fe(m+1)*n) ((-1)^m/W(m+1)) ];
    end
    % rownanie wspol c
    c = A\H;
    h = c(1:L); sigma=abs(c(M));
    g = h'/2; g(1)=2*g(1); g = [fliplr(g(2:L)) g];
    figure(59);
    sgtitle('Odp impulsowa i czestoliwosciowa filtra');
    subplot(221); stem(h); title('Odp impulsowa - polowa');
    subplot(222); stem(g); title('Odp impulsowa - cala');
    
    for k=0:K-1
        H(k+1) = cos(pi*fz(k+1)*n)* h;
        Herr(k+1) = Wg(k+1) * (H(k+1) - Hd(k+1));
    end
    subplot(223); plot(fz, Hd, 'r', fz, H, 'b'); grid;
    title('Charakterystyka czestotliwosciowa'); ylabel('Ampl.');
    subplot(224); plot(fz, Herr); grid;
    title('Charakterystyka czestotliwosciowa - blad'); ylabel('Ampl.');
    
    % szukanie ekstremow
    Hmax= []; i_maximum = [];
    for p = 1 : 2
        if (p==1) k1=k11; k2=k12; end
        if (p==2) k1=k21; k2=k22; end
        Hmax = [Hmax Herr(k1)]; i_maximum = [i_maximum k1];
        k = k1 + 1;
        while(Herr(k-1) == Herr(k)) k = k + 1; end
        if (Herr(k) < Herr(k+1))
            sgn = 1;
        else
            sgn = -1;
        end
        k=k+1;
        while(k<=k2)
            if (sgn==1)
                while( (k<k2) & (Herr(k-1) < Herr(k))) k = k + 1; end
            end
            if (sgn==-1)
                while( (k<k2) & (Herr(k-1) > Herr(k))) k = k + 1; end
            end
            sgn = -sgn;
            Hmax = [Hmax Herr(k)]; i_maximum = [i_maximum k];
            k=k+1;
        end
    end
    figure(60);
    subplot(211);
    plot(fz(i_maximum),Hmax,'or',fz,Herr,'b');
    grid; title('B³¹d charakterystyki i jego ekstrema');
    
    % Wybranie M+1 najwiekszych
    if length(Hmax)>M
        IM = []; G = abs(Hmax); LenG = length(G);
        while( LenG > 0)
            Gmx = max(G); imx = find(G==Gmx);
            LenGmx = length(imx);
            IM = [ IM i_maximum(imx)];
            G(imx) = 0; LenG = LenG - LenGmx;
        end
        IM = IM(1:M); IM = sort(IM); i_maximum = IM;
    end
    sigmaMAX = max(abs(Hmax));
    feMAX = fz(i_maximum);
    subplot(212);
    plot(fz(i_maximum),Herr(i_maximum),'or',fz,Herr,'b');
    grid; title('B³¹d charakterystyki i M+1 najwiêkszych ekstremów'); pause(5);
    
end

fz = fz/2;
figure(61);
subplot(211); stem(g); title('Wynikowa odp impulsowa filtra');
subplot(212); plot(fz(i_maximum),Herr(i_maximum),'or',fz,Herr,'b'); grid;
title('B³¹d H(f) + jego EKSTREMA');
figure(62);
subplot(211);
plot(fz,Hd,'r',fz,H,'b'); grid; title('Wynikowe H(f)');
subplot(212);
plot(fz,20*log10(H),'b'); grid; title('Wynikowe H(f) w dB');

%% Æwiczenie: Projektowanie nierekursywnych filtrów cyfrowych metod¹
% okien z zastosowaniem okna Kaisera
% Podaj parametry filtra (np. pasmowozaporowego)
 fpr = 1000;	% czêstotliwoœæ próbkowania [Hz]
 fd1 = 150;	% czêstotliwoœæ dolna 1 [Hz]
 fd2 = 200;	% czêstotliwoœæ dolna 2 [Hz]
 fg1 = 300;	% czêstotliwoœæ górna 1 [Hz]
 fg2 = 350;	% czêstotliwoœæ górna 2 [Hz]
 dp = 0.001;	% oscylacje w paœmie przepustowym np. 0.1, 0.01, 0.001
 ds = 0.0001;	% oscylacje w paœmie zaporowym np. 0.001, 0.001, 0.0001
 typ = ["LowPass", "HighPass", "BandPass","BandStop"];
 % lp=LowPass, hp=HighPass, bp=BandPass, bs=BandStop
 
 for i=1:length(typ)
    if i==1
        df=fd2-fd1;     % filtr low pass
        fc=((fd1+fd2)/2)/fpr;
        wc=2*pi*fc;
    end
    if i==2
        df=fg2-fg1;     % filtr high pass
        fc=((fg1+fg2)/2)/fpr;
        wc=2*pi*fc;
    end
    if (i==3 || i==4)
        df1=fd2-fd1; df2=fg2-fg1;   % filtry pasmowe
        df = min(df1, df2);
        f1 = (fd1+df/2)/fpr;
        f2 = (fg2-df/2)/fpr;
        w1 = 2*pi*f1; w2 = 2*pi*f2;
    end
    d = min(dp, ds); A = -20*log10(d);
    % wzory na okno Kaisera w zale¿nosci od tlumienia
    if (A>=50) beta = 0.1102*(A-8.7);end
    if (A>21 & A<50) beta = (0.5842*(A-21)^0.4)+0.07886*(A-21);end
    if (A<=21) beta = 0; end
    if (A>21) D = (A-7.95)/14.36; end
    if (A<=21) D = 0.922; end
    % Wyznaczenie dlugosci
    N = ceil((D*fpr/df)+1); if (rem(N,2)==0) N=N+1; end
    M = (N-1)/2; m =1:M; n=1:N; % wektory indeksowe
    % Generacja okna czaswego
    temp =  beta * sqrt(1-((n-1-M).^2./M.^2));
    wb =  besseli( 0, temp ) / besseli(0,beta);
    figure(63+2*i-1); subplot(311); plot(n, wb); grid; 
    title(["Okno czas. Kaisera M=" + num2str(M) " Filtr: " + typ(i)]);
    
    % Tworzenie odp. impulsowej, wzory z tabelki
    if (i==1) h=2*fc*sin(wc*m)./(wc*m); h=[ fliplr(h) 2*fc h]; end % filtr LP
    if (i==2) h=-2*fc*sin(wc*m)./(wc*m); h=[ fliplr(h) 1-2*fc h]; end % filtr HP
    if (i==3) % filtr BP
        h = 2*f2*sin(w2*m)./(w2*m) - 2*f1*sin(w1*m)./(w1*m);
        h = [ fliplr(h) 2*(f2-f1) h];
    end
    if (i==4) % filtr BS
        h = 2*f1*sin(w1*m)./(w1*m) - 2*f2*sin(w2*m)./(w2*m);
        h = [ fliplr(h) 1+2*(f1-f2) h];
    end
    subplot(312); plot(n,h); grid;
    title(["Odp impulsowa filtra: " + typ(i)]);
    
    % Wymono¿enie odp impulswoej z oknem
    hw = h.*wb;
    subplot(313); plot(n,hw,'b'); grid;
    title(["Iloczyn okna i odp. impulsowej filtra: " + typ(i)]);
    
    % Charakterystyka czestotliwosciowa
    NF = 1000; fmin = 0; fmax = fpr/2;      % wartoœci parametrów charakterystyki
    f = fmin : (fmax-fmin)/(NF-1) : fmax;	% czêstotliwoœæ
    w = 2*pi*f/fpr;                         % pulsacja
    HW = freqz(hw,1,w);
    figure(63+2*i);
    sgtitle(['Charakterystyka czestotliwosciowa filtra: ' typ(i)]);
    subplot(211); plot(f, abs(HW)); grid; ylabel('Ampl.');
    xlabel('Czestotliwoœæ [Hz]');
    title('Modu³ odp. czêstotliwoœciowej');
    subplot(212); plot(f, unwrap(angle(HW))); grid; ylabel('rad.');
    xlabel('Czestotliwoœæ [Hz]');
    title('Faza odp. czêstotliwoœciowej');

 end

 %% Algorytm interpolacji sygna³u za pomoc¹ dyskretnej transformacji Fouriera DFT (FFT)
 M=24; N=16; n=0:N-1; x=sin(2*pi/8*n);
 X = fft(x);    % obliczenie fft sygnalu rzeczywistego
 % Wstawienie do widma w srodku symetrii zer, skrajnie zer obliczyæ
 % jak 0.5 bo sk³adowa zosta³¹ "rozbita"
 X = [ X(1:N/2) 0.5*X(N/2+1) zeros(1,M-N-1) conj(0.5*X(N/2+1)) X(N/2+2:N)];
 y = M/N*real(ifft(X));
 figure(72); sgtitle('Filtr interpolujacy');
 subplot(211); stem(x); title('Orginalny'); 	% sygna³ wejœciowy
 subplot(212); stem(y); title('Interpolowany'); % interpolowany sygna³ wejœciowy

 %% Æwiczenie: Projektowanie specjalnych filtrów cyfrowych metod¹ okien
 % Parametry
 M = 20; N = 2*M+1; n=1:M;
 typ = ["Hilberta", "ró¿niczkuj¹cy", "interpoluj¹cy"];
 
 for i=1:3
     % Generowanie polowy odpowiedzi impulsowej
    if(i==1)
        h = 2/pi*sin(pi*n/2).^2 ./ n;   % Odp impulsowa filtra Hilberta
    end
    if i==2
        h = cos(pi*n)./n;               % Odp impuslowa filtra rozn.
    end
    if i==3
        K = 5; wc=pi/K; fc=wc/(2*pi);
        h = 2*fc*sin(wc*n)./(wc*n);     % ODp impuls filtra interp.
    end
    if (i==1 || i==2)
        h = [-h(M:-1:1) 0 h(1:M)];      % Odbicie odpowiedzi
    else
       h = K*[-h(M:-1:1) 2*fc h(1:M)];  % Odbicie i skalowanie przez K
    end
    
    % Mno¿enie odp z oknem
    w = blackman(N)';
    hw = h.*w;
    
    % Widmo Fouriera
    m = -M:1:M;
    NF = 500; fn=0.5*(1:NF-1)/NF;
    for k=1:NF-1
        H(k)=sum( h .* exp(-j*2*pi*fn(k)*m) );
        HW(k)=sum( hw .* exp(-j*2*pi*fn(k)*m) );
    end
    
    % wykresy
    figure(72+2*i-1);
    sgtitle(["Odpowiedz impulsowa filtra: " + typ(i)]);
    subplot(211); stem(m,h); grid; title('h(n)'); xlabel('n');
	subplot(212); stem(m,hw); grid; title('hw(n) - wymnozenie z oknem');
    xlabel('n');
    
    figure(72+2*i);
    % Zastosowanie filtrów 1 i 2
    if(i<3)
       % Sygna³ testowy
       Nx=200; fx=50; fpr=1000; n=0:Nx-1; x=cos(2*pi*fx/fpr*n);
       y = conv(x, hw); % filtracja sygna³u odpowiedzi¹
       % odciêcie stanów przejœciowych (po N?1 próbek) z przodu i z ty³u sygna³u yz(n)
       yp = y(N:Nx);
       % odciêcie tych próbek w xz(n), dla których nie ma poprawnych odpowiedników w yz(n)
       xp = x(M+1:Nx-M);
       if (i==1)
           z = xp + 1i*yp;  % sygnal analityczny
           Ny = ceil(fpr/fx); k=1:Ny;
           subplot(311); plot(k, xp(k), 'b', k, yp(k), 'r');
           title('Sygnal analityczny'); legend('real(z)', 'imag(z)');
           subplot(312); plot(xp, yp); title('Imag(Real(z))'); grid;
           subplot(313); plot(abs(fft(z)));
           title('Widmo sygna³u analitycznego'); grid;
       else
           % Filtr rozniczkujacy
           Ny = ceil(fpr/fx); k=1:Ny;
           plot(k, xp(k), 'b', k, yp(k), 'r');
           title('Sygnal filtrowany filtrem rozczkujacym');
       end
    end
    if (i==3)
         % generacja sygna³u testowego x(n)
        Nx=50; fx=50; fpr=1000; n=0:Nx-1; x=cos(2*pi*fx/fpr*n);
        xz=[]; KNx=K*Nx; xz=zeros(1,KNx); xz(1:K:KNx)=x(1:Nx); % dodanie zer
        % filtracja xz(n) za pomoc¹ odp. impulsowej hw(n); otrzymujemy Nx+N?1 próbek
        yz=conv(xz,hw);
        % odciêcie stanów przejœciowych (po N?1 próbek) z przodu i z ty³u sygna³u yz(n)
        yp=yz(N:KNx);
        % odciêcie tych próbek w xz(n), dla których nie ma poprawnych odpowiedników w yz(n)
        xp=xz(M+1:KNx-M);
        Ny=length(yp); k=1:Ny;
        plot(k,xp(k),'or',k,yp(k),'-b');
        title('Sygna³ filtrowany filtrem interpolujacym'); grid; % porównanie
    end

 end
 
