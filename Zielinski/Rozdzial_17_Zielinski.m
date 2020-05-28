% Rodzial 17 Zielniski
%                       Mateusz Krupnik

% Æwiczenie: Czasowo-czêsotliwoœciowa reprezentacja Gabora
clc; clear all; close all;
% Parametry wejœciowe
nw = 64;    % D³ugoœæ okna
L = 128;    % Dlugoœæ okna po uzupe³nieniu zerami
dM = 4;     % Przesuniecie w czasie
dN = 4;     % Przesuniecie w czestotliwoœci
% dMdN < L !
% Okno czasowe do analizy (musi t³umiæ sygna³ jak np. Gaussa)
w = blackman(nw)';  % Okno Blackmana
% Uzuep³nienie zerami po bokach do d³ugoœci L
w = [zeros(1, (L-nw)/2) w zeros(1, (L-nw)/2) ];

% Genreacja odwrotengo okna analizy
M = L/dM; N = L/dN; ww = [ w w ];
H = []; k = 0 : L-1;
for p=0:dM-1
   for q=0:dN-1
       h = ww(k+q*N+1) .* exp(-1i*2*pi*p*k/dM); % Nowa funkcja bazy
       H = [H; h];      % Dodanie nowej funkcji bazy do bazy
   end
end

mi = zeros(dM*dN, 1); mi(1,1)=dM/N;
dw = pinv(H)*mi; dw=dw';

figure(1)
sgtitle(['Czasowo-czêsotliwoœciowa reprezentacja Gabora']);
subplot(211); plot(real(w)); title('Okno czasowe Blackmana z zerami'); grid;
subplot(212); plot(real(dw)); title('Okno dualne do powy¿szego'); grid;

% Generacja sygna³u zmiennego w czestotliwosci z czasem
% Parametry
fpr = 128; f0 = 0; df = 24; dt = 1/fpr; n = 0:L-1; t = n*dt;
x = sin(2*pi*(f0*t+0.5*df*t.^2));

% Analiza
% Punkty czasu L/dM, punkty czest. L/dN
dwz = [dw zeros(1,2*L)];
for k=0:2*L/dM-1
   okno = dwz(L+1:L+L);             % Wyciecie okna
   widmo = fft(x.*okno);            % FFT
   tf(k+1, :) = widmo(1:dN:L);      % Decymacja co krok czest.
   dwz = [zeros(1,dM) dwz(1:3*L-dM)];   % przesuniecie okna o krok czas.
end
figure(2)
sgtitle(['Czasowo-czêsotliwoœciowa reprezentacja Gabora']);
subplot(2,2,[1, 2]); plot(t, x); title('Sygna³ wejœciowy'); xlabel('Czas [s]');
subplot(223); mesh( abs(tf')); title('Wykres 3D'); grid;
subplot(224); contour( abs(tf')); title('Kontur'); grid;

% Synteza sygnalu
wz = [w zeros(1,2*L)];  % Ono syntezy uzupelnione zerami
y = zeros(1, L);
temp = (2*pi/L)*dN;
n = 0:1:L-1;
for m=0:1:2*L/dM-1
   ww = wz(L+1:L+L);    % Okno z zerami
   for k=0:1:L/dN-1
      y = y + tf(m+1, k+1)*(ww .* ( cos(temp*k*n) + 1i*sin(temp*k*n))); 
   end
   wz = [ zeros(1, dM) wz(1:3*L-dM)];   % Przesuniecie okna o krok czas.
end
y = real(y);

figure(3); sgtitle(['Czasowo-czêsotliwoœciowa reprezentacja Gabora']);
subplot(311); plot(n, y); title('Sygna³ wyjœciowy z syntezy'); grid;
subplot(312); plot(n, x, n, y); title('Porownaie z wejœciowym'); grid;
legend('WEJ x(t)','WYJ y(t)');
subplot(313); plot(n, y-x); title('B³ad y-x'); grid;

blad_max = max(abs(y-x))
blad_std = std(x-y)

%% Æwiczenie: Krótkoczasowa transformacja Fouriera
% Parametry wejœciowe
M=32;       % po³owa d³ugoœci okna (ca³e okno N=2M-1)
Nx=128;     % d³ugoœæ sygna³u testowego
% Sygna³ testowy z modulacj¹ czêstotliwoœci typu LFM i SFM
fpr=128; f0=0; df=32; fn=16; fm=3; dfm=12; dt=1/fpr; n=0:Nx-1; t=n*dt;
x1 = sin(2*pi*(f0*t+0.5*df*t.^2));
x2 = sin( 2*pi* (fn*t + (dfm/(2*pi*fm))*sin(2*pi*fm*t)) );
% Analiza TF ? krótkoczasowa reprezentacja Fouriera
x1 = hilbert(x1); x2 = hilbert(x2);     % Sygna³ analityczny 
w = hanning(2*M-1)';                    % Okno Hanninga
for n = M : Nx-M+1  
    xx1 = x1(n-(M-1): 1 :n+(M-1)); xx1 = xx1 .* w; xx1 = [ xx1 0 ];
    X1(:,n-M+1) = fftshift(abs(fft(xx1))');
    xx2 = x2(n-(M-1): 1 :n+(M-1)); xx2 = xx2 .* w; xx2 = [ xx2 0 ];
    X2(:,n-M+1) = fftshift(abs(fft(xx2))');
end
% Rysunek widma TF
t_=t(M:Nx-M+1); f=fpr/(2*M)*(-M:M-1);
figure(4); sgtitle('Krótkoczasowa transformacja Fouriera sygna³u LFM');
subplot(2,2, [1,2]); plot(t, x1); title('Synga³ LFM'); grid;
subplot(223); mesh(t_,f,X1); view(-40,70); axis tight;
title('Wykres analizy w 3D');
xlabel('Czas [s]'); ylabel('Czêstotliwoœæ [Hz]');
subplot(224); imagesc(t_,f,X1); title('Rzut z góry');
xlabel('Czas [s]'); ylabel('Czêstotliwoœæ [Hz]');

figure(5); sgtitle('Krótkoczasowa transformacja Fouriera sygna³u SFM');
subplot(2,2, [1,2]); plot(t, x2); title('Synga³ SFM'); grid;
subplot(223); mesh(t_,f,X2); view(-40,70); axis tight;
title('Wykres analizy w 3D');
xlabel('Czas [s]'); ylabel('Czêstotliwoœæ [Hz]');
subplot(224); imagesc(t_,f,X2); title('Rzut z góry');
xlabel('Czas [s]'); ylabel('Czêstotliwoœæ [Hz]');

%% Æwiczenie: Generacja funkcji skaluj¹cych i falek.
clear all;
niter = 10; % liczba iteracji
c = 0; d = 1;      % {c=1, d=0} ? funkcja skaluj¹ca, {c=0, d=1} ? falka
% definicja wspó³czynników filtrów h0 i h1 systemu falkowego Db4 (17.62) (17.57)
h0 = [ (1+sqrt(3))/(4*sqrt(2)) (3+sqrt(3))/(4*sqrt(2)) ...
    (3-sqrt(3))/(4*sqrt(2)) (1-sqrt(3))/(4*sqrt(2)) ];
N = length(h0); n = 0:N-1;
h1 = (-1).^n .* h0(N:-1:1);
% synteza ? wed³ug schematu drzewa filtrów z rysunku 17.15
c = [ 0 c 0 ]; % aproksymacje ?0
d = [ 0 d 0 ]; % detale ?0
c = conv(c,h0) + conv(d,h1);
for n = 1 : niter
    for k = 1:length(c)
        c0(2*k-1) = c(k);
        c0(2*k) = 0;
    end
    c0 = [ 0 c0 ];
    c = conv(c0, h0);
end
figure(6);
plot(c); title('Przyk³adowa falka');

%% Æwiczenie: Transformacja falkowa
clear all;
% Parametry programu
niter = 3;          % liczba iteracji
nx = 2^niter*32;	% d³ugoœæ sygna³u
% Definicja wspó³czynników filtra LP syntezy h0s, np. Db4
h0s = [ (1+sqrt(3))/(4*sqrt(2)) (3+sqrt(3))/(4*sqrt(2)) ...
    (3-sqrt(3))/(4*sqrt(2)) (1-sqrt(3))/(4*sqrt(2)) ];
% Oblicz pozosta³e filtry
N = length(h0s); n = 0:N-1;
h1s = (-1).^n .* h0s(N:-1:1);       % filtr HP syntezy
h0a = h0s(N:-1:1); h1a=h1s(N:-1:1); % filtry LP i HP analizy
% Sygna³ testowy
x1 = sin(2*pi*(1:nx)/32);
x2 = rand(1,nx);
% Analiza
cc1 = x1;
cc2 = x2;
for m=1:niter
    c01 = conv(cc1,h0a);        % filtracja LP x1
    d01 = conv(cc1,h1a);        % filtracja HP x1
    
    k1=N:2:length(d01)-(N-1); kp1=1:length(k1);
    ord1(m)=length(kp1); dd1(m,kp1) = d01( k1 );
    k1=N:2:length(c01)-(N-1); cc1=c01( k1 );
    
    c02 = conv(cc2,h0a);        % filtracja LP x2
    d02 = conv(cc2,h1a);        % filtracja HP x2
    
    k2=N:2:length(d02)-(N-1); kp2=1:length(k2);
    ord2(m)=length(kp2); dd2(m,kp2) = d02( k2 );
    k2=N:2:length(c02)-(N-1); cc2=c02( k2 );
end
% Synteza sygna³ów
c1=cc1; c2=cc2;
for m=niter:-1:1
    c01=[]; d01=[];
    for k = 1:length(c1)
        c01(2*k-1)=c1(k); c01(2*k)=0;
    end
    c1 = conv(c01,h0s); nc1=length(c1);
    for k = 1:ord1(m)
        d01(2*k-1) = dd1(m,k); d01(2*k) = 0;
    end
    d1 = conv(d01,h1s); nd1=length(d1);
    c1 = c1(1:nd1);
    c1 = c1 + d1;
    
    c02=[]; d02=[];
    for k = 1:length(c2)
        c02(2*k-1)=c2(k); c02(2*k)=0;
    end
    c2 = conv(c02,h0s); nc2=length(c2);
    for k = 1:ord2(m)
        d02(2*k-1) = dd2(m,k); d02(2*k) = 0;
    end
    d2 = conv(d02,h1s); nd2=length(d2);
    c2 = c2(1:nd2);
    c2 = c2 + d2;
end

% Wykresy koñcowe
n1 = 2*(N-1)*niter : length(c1)-2*(N-1)*niter+1;
figure(7); sgtitle('Analiza sygna³u sinusoidalnego');
subplot(311); title('Sygna³ wejsciowy');
plot(x1); grid;
subplot(312);
plot(n1,c1(n1)); title('Sygna³ wyjœciowy'); grid;
subplot(313); title('Blad analizy'); grid; plot(n1, x1(n1)-c1(n1));

n2 = 2*(N-1)*niter : length(c2)-2*(N-1)*niter+1;
figure(8); sgtitle('Analiza sygna³u losowego');
subplot(311); title('Sygna³ wejsciowy');
plot(x2); grid;
subplot(312);
plot(n2,c2(n2)); title('Sygna³ wyjœciowy'); grid;
subplot(313); title('Blad analizy'); grid; plot(n2, x2(n2)-c2(n2));
