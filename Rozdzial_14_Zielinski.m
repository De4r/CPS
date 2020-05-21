% Rodzial 13 Zielniski
%                       Mateusz Krupnik

% Æwiczenie: Filtry adaptacyjne typu LMS (NLMS) ? losowego gradientu
% Oznaczenia: x - syg filtrowany, d - sygnal odniesienia
%             y - przefiltrowany adaptacyjnie sygnal x,
%             e = d - y - sygnal bledu adaptacji
clear all; close all; clc;
% Wybor okna
okno_ = 2; % 1 - brak, 2-gaussa, 3-alfa*t, 4-exp(-alfa*t)
alg =2;     %1 - LMS, 2 - NLMS
algo = ["LMS", "NLMS"];
% Parametry filtrów
% LMS
M = 50; mi = 0.1;   % mi <1 i >0
% NLMS
eng = 0.0; beta = 1 - 1/M;  % energia poczatkowa, wspolczynnik pamieci
gamma = 0.001;              % stala bezpieczenstwa mianownika!

% GEneracja sygnalu testowego
Nx = 1000;  % probki
fpr = 1000; % Hz
A = 1;      % Ampl
f0 = 0;     % czestotliwosc poczatkowa sygnalu
df = 25;    % przyrost czest. na sek
dt = 1/fpr; t=0:dt:(Nx-1)*dt;
% Sygna³
s = A*cos(2*pi*(f0*t + df/2*t.^2));
% okna obwiedniowe
if (okno_==2) alfa=10; w=exp(-alfa*pi*(t-0.5).^2);end
if (okno_==3) alfa=5; w=alfa*t; end
if (okno_==4) alfa=5; w=exp(-alfa*t); end
if (okno_~=1) s = s.* w; end
okno = ["Brak", "Gaussa", "Liniowe", "Wyk³adnicze"];

for i=1:2
    if i==1
        % KAsowanie interferencji sieci - sinusida 50Hz przesunieta w fazie
        P = 0; % brak predykcji
        x = 0.1*sin(2*pi*50*t-pi/5); % sieæ przesuniêta w fazie
        d = s + 0.5*sin(2*pi*50*t); % sygna³ + sieæ
    else
        % Odszumianie sygnalu z szumu normlnaego
        P = 1; % rz¹d predykcji (do zmiany: 1,2,3,...)
        x = s + 0.25*randn(1,Nx); % sygna³ + szum
        d = [ x(1+P:length(x)) zeros(1, P) ]; % odniesieniem sygna³ "przyspieszony" o P próbek
    end
    tempp = " przyk³ad: " + num2str(i) + " Okno: " + okno(okno_) + " Alg: " + algo(alg);
    figure();
    sgtitle("Sygna³y wejœciowe" + tempp);
    subplot(211); plot(t, x); grid;
    title('Syg. wej.: x(t)'); xlabel('Czas [s]'); ylabel('Ampl.');
    subplot(212); plot(t, d); grid;
    title('Syg. wej.: d(t)=s(t)+x(t)'); xlabel('Czas [s]'); ylabel('Ampl.');
    
    % Filtracja adaptacyjna, 4 wiersze dla 4 sygnalow
    bx = zeros(1, M);   % buffor wejsciowy x
    h = zeros(1, M);    % Wagi filtra LMS
    y = []; e = [];     % wektory wyniku
    for k=1:length(x)
        bx = [x(k) bx(1:M-1)];
        dest = h * bx';

        err = d(k) - dest;
        if (alg==1)
            h = h + (2*mi*err*bx);
        else
            eng = bx*bx';
            h = h + ( (2*mi)/(gamma + eng) * err * bx);
        end
        y = [y dest];
        e = [e err];
    end
    % Wykresy wyników
    figure();
    sgtitle("Sygna³y wyj.," + tempp);
    subplot(211); plot(t, y); grid;
    title('Syg. wyj.: y(t)'); xlabel('Czas [s]'); ylabel('Ampl.');
    subplot(212); plot(t, e); grid;
    title('Syg. wyj.: e(t)=d-y'); xlabel('Czas [s]'); ylabel('Ampl.');
   
    figure()
    if (i==1)
        subplot(111); plot(t,s,'r',t,e,'b'); grid; xlabel('czas [sek]');
    end
    if (i==2)
        n=1:Nx-P;
        subplot(111); plot(t(n),s(n),'r',t(n),y(n),'b'); grid; xlabel('czas [sek]');
    end
    title("Orgina³ (czerwony) i wynik filtracji (niebieski)" + tempp);
end