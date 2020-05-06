% Lab 6. Dla zadanych transmitancji wykreœliæ ich
% odpowiedzi impulsowe, skokowe, charakterystyki 
% amplitudowe i fazowe oraz zbadac stabilnoœæ.
% Dla filtra FIR sprawdzic jego dzia³anie.
%                                   Mateusz Krupnik

%%%%%%%%%%%%% CZESC 1.
clc; clear all; close all;
L=[0.0675 0.1349 0.0675; 0.0412 0.0824 0.0412;...
    0.0996 0.1297 0.0996; 0.1239 0.0662 0.1239]; % Licznik
M=[1 -1.1430 0.4128; 1 -1.4409 0.6737; 1 -1.6099 0.6794; ...
    1 -1.4412 0.6979]; % Mianownik
d_k = zeros(1, 1001); d_k(1,1)=1;
fs=1000; f=100; t=0:(1/fs):1;
for i=1:size(L,1)
    disp(['System nr: ' num2str(i)]);
    l = (L(i, :)); m = (M(i, :));
    sys=tf(l, m) % Transfer Function
    % Wykresy systemu
    [h, w] = freqz(l, m);	% Char. ampl i fazowa
    Mag = 10*log(abs(h));     % Amplituda w skali log
    F = phase(h)*180/pi;	% Faza w stopiach
    w = w/pi;               % skalowanie
    figure(4*i-3)
    sgtitle(['Charakterystyki dla systemu nr: ' num2str(i)]);
    subplot(211); plot(w, Mag); 
    ylabel('Amplituda [dB]'); xlabel('Czest. znorm.');
    subplot(212); plot(w, F); 
    ylabel('Faza [\circ]'); xlabel('Czest. znorm.');
    figure(4*i-2)
    dimpulse(l, m);         % impuls dla Z transmitacji
    
    y = filter(l, m, d_k);    % odpowiedz filtra o Z transmitacji na x
    figure(4*i-1)
    plot(t,d_k, t, y);
    title(['OdpowiedŸ dla systemu nr: ' num2str(i)]); xlabel('Czas [s]');
    ylabel('Amplituda'); legend('Wymuszenie', 'OdpowiedŸ');
    
    a=0; b=0; r=1; % 
    x = linspace(a-r,a+r,100);
    y1=sqrt(r^2-(x-a) .^2)+b;
    y2=-sqrt(r^2-(x-a) .^2)+b;
    figure(4*i)
    plot(x,[y1; y2],'b'); grid on; axis equal; hold on;
    pzmap(l, m); % Wyrysowanie zer i biegunow
    [z, p, k] = tf2zpk(l, m);
    
end

%% Sprawdzenie dzia³ania filtrów cyfrowych
% Dane sygna³ów
f1=100; f2=250; f3=400; fs=1000;    % Czestotliwosci skladowych i Nyquista
A1=1; A2=0.8; A3=0.65;              % Aplitudy skladowych
t=0:(1/fs):1.023;                   % Wektor czasu
% Sygna³ wymuszenia - sk³adowa 3 harmonicznych
x=A1*sin(2*pi*f1*t)+A2*sin(2*pi*f2*t)+A3*sin(2*pi*f3*t);

%% Filtr dolnoprzepustowy
% Parametry
fo1 = 300; fo2 = 200;           % Czestotliwosci ociecia
wn1 = fo1*2/fs; wn2 = fo2*2/fs; % Znormalizowane czestotliwosci
wn = [wn1, wn2]; fo = [fo1, fo2];
% Filtr Butter - dolnoprzepustowy
N = [2, 4 ,8];
for i=1:length(wn)
    [l, m] = butter(N(1,3), wn(i));     % Rz¹d oraz czest. odciecia
    y = filter(l, m, x);                % Odpowiedz filtra
    % Obliczenie transformaty Fouriera 
    [f_w, Moc, Wid] = fft_from_signal([x; y], fs);   % Funkcja z Lab 4 i 5
    % Dzwiek
    sound(x); pause(t(end)); sound(y); pause(t(end));
    % Wykresy
    figure(16+i); sgtitle(['Filtr dolnoprzepustowy f_o=' num2str(fo(i))]);
    subplot(321); plot(t, x); xlabel('Czas [s]'); ylabel('Amplituda'); grid;
    title('Sygna³ wymuszaj¹cy');
    subplot(322); plot(t, y); xlabel('Czas [s]'); ylabel('Amplituda'); grid;
    title('Sygna³ po filtracji');
    subplot(323); plot(f_w, Wid(1,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
    title('Widmo wymuszenia');
    subplot(324); plot(f_w, Wid(2,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
    title('Widmo odpowiedzi');
    subplot(325); plot(f_w, Moc(1,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
    title('Moc widmowa wymuszenia');
    subplot(326); plot(f_w, Moc(2,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
    title('Moc widmowa odpowiedzi');
end

%% Filtr srodkowoprzepustowy
% Parametry z poprzedniego filtru
[l, m] = butter(N(1,3), fliplr(wn));    % Odwrócenie kolejnoœci cz. odc.
y = filter(l, m, x);                    % Odpowiedz filtra
[f_w, Moc, Wid] = fft_from_signal([x; y], fs);   % Funkcja z Lab 4 i 5

% Dzwiek
sound(x); pause(t(end)); sound(y); pause(t(end));
% Wykresy
figure(19); sgtitle(['Filtr srodkowoprzepustowy f_o1=' num2str(fo(2)) ' f_o2=' num2str(fo(1))]);
subplot(321); plot(t, x); xlabel('Czas [s]'); ylabel('Amplituda'); grid;
title('Sygna³ wymuszaj¹cy');
subplot(322); plot(t, y); xlabel('Czas [s]'); ylabel('Amplituda'); grid;
title('Sygna³ po filtracji');
subplot(323); plot(f_w, Wid(1,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Widmo wymuszenia');
subplot(324); plot(f_w, Wid(2,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Widmo odpowiedzi');
subplot(325); plot(f_w, Moc(1,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Moc widmowa wymuszenia');
subplot(326); plot(f_w, Moc(2,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Moc widmowa odpowiedzi');

%% Filtr górnoprzepustowy
% Parametry z poprzednich filtrów
[l, m] = butter(N(1,3), wn(1), 'high');
y = filter(l, m, x);                    % Odpowiedz filtra
[f_w, Moc, Wid] = fft_from_signal([x; y], fs);   % Funkcja z Lab 4 i 5

% Dzwiek
sound(x); pause(t(end)); sound(y); pause(t(end));
% Wykresy
figure(20); sgtitle(['Filtr górnoprzepustowy f_o1=' num2str(fo(1))]);
subplot(321); plot(t, x); xlabel('Czas [s]'); ylabel('Amplituda'); grid;
title('Sygna³ wymuszaj¹cy');
subplot(322); plot(t, y); xlabel('Czas [s]'); ylabel('Amplituda'); grid;
title('Sygna³ po filtracji');
subplot(323); plot(f_w, Wid(1,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Widmo wymuszenia');
subplot(324); plot(f_w, Wid(2,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Widmo odpowiedzi');
subplot(325); plot(f_w, Moc(1,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Moc widmowa wymuszenia');
subplot(326); plot(f_w, Moc(2,:)); xlabel('Czestotliwoœæ [Hz]'); ylabel('Amplituda'); grid;
title('Moc widmowa odpowiedzi');

%% Projekt filtra i jego dzialanie
% Parametry, niektóre z poprzednich filtrów
format long e;  % Typ formatowania zmiennych, 16 miejsc, wykladniczo
F_=[0 0.1 0.2 0.5 0.7 1];
M_=[1 1 1 0 0 0];

[l, m] = yulewalk(N(1,3), F_, M_);      % Metoda Yule-Walkera
b = fir2(128, F_, M_); fir2(128, F_, M_);


[h, w] = freqz(b, 1); Mag = 20*log(abs(h)); Fi = phase(h)*180/pi; w = w/pi;
[h, w_] = freqz(l, m); Mag_ = 20*log(abs(h)); Fi_ = phase(h)*180/pi; w_ = w_/pi;
% Wykresy
figure(21);
sgtitle('Porownanie roznych metod projektowania');
subplot(211); plot(w, Mag, w_, Mag_); legend('FIR2', 'Yule-Walker'); grid;
xlabel('Czest. znorm.'); ylabel('Amplituda [dB]');
title('Charakterystyka amplitudowa');
subplot(212); plot(w, Fi, w_, Fi_); legend('FIR2', 'Yule-Walker'); grid;
xlabel('Czest. znorm.'); ylabel('Faza [\circ]');
title('Charakterystyka fazowa');

