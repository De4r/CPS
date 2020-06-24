% Lab 5. Analiza uk³adów RC i RLC
%                                   Mateusz Krupnik
%% Analiza uk³adów RC
% Dane uk³adu
R=1000; C=10^(-6);
L=[1]; M=[(R*C) 1]; % Licznik, Mianownik
sys=tf(L,M) % Transfer function (licznik, mianownik)
% Wykresy analizy uk³adu
figure(1)
freqs(L,M)      % analiza amplitudowo i fazowo - czestotliwosciowa
figure(2)
impulse(L,M)	% odpowiedz impulsowa ukladu
figure(3)
step(L,M)       % odpowiedz skokowa ukladu
figure(4)
iopzplot(sys)   % wykres zer i biegunow dla ukladow wej/ wyj
[z,p,k]=tf2zp(L,M)  % konwersja Transfer Function na zera i bieguny

%% Analiza uk³adu RLC
% Dane ukladu
R=1000; C=10^(-6); Li=1;
L=[1]; M=[(Li*C) (R*C) 1]; % licznik i mianownik
sys=tf(L,M)	% Transfer function (licznik, mianownik)
% Wykresy analizy uk³adu
figure(5)
freqs(L,M)      % analiza amplitudowo i fazowo - czestotliwosciowa
figure(6)
impulse(L,M)    % odpowiedz impulsowa ukladu
figure(7)
step(L,M)       % odpowiedz skokowa ukladu
figure(8)
iopzplot(sys)   % wykres zer i biegunow dla ukladow wej/ wyj
[z,p,k]=tf2zp(L,M)  % konwersja Transfer Function na zera i bieguny