% Lab 6. Dla zadanych transmitancji wykreœliæ ich
% odpowiedzi impulsowe, skokowe, charakterystyki 
% amplitudowe i fazowe oraz zbadac stabilnoœæ.
%                                   Mateusz Krupnik

% CZESC 1: Wykreœlenie charakterystyk dla 4 zestawów wspó³czynnikow
% transmitacji uk³adów (KODY 1-4 w Lab. 6)
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
    [h, w] = freqz(l, m, fs);	% Char. ampl i fazowa
    Mag = 20*log10(abs(h));   % Amplituda w skali log
    F = phase(h)*180/pi;	% Faza w stopiach
    w = w/pi;               % skalowanie
    figure(4*i-3)
    sgtitle(['Charakterystyki dla systemu nr: ' num2str(i)]);
    subplot(211); plot(w, Mag); 
    ylabel('Amplituda [dB]'); grid;
    xlabel('Czestoœæ znormalizowa');
    subplot(212); plot(w, F); 
    ylabel('Faza [\circ]'); grid;
    xlabel('Czêstoœæ znormalizowa');
    figure(4*i-2)
    subplot(211);
    dimpulse(l, m); grid;	% impuls dla Z transmitacji
    subplot(212);
    dstep(l, m); grid;      % odp skokowa
    
    y = filter(l, m, d_k);    % odpowiedz filtra o Z transmitacji na x
    figure(4*i-1)
    plot(t(1:50),d_k(1:50), t(1:50), y(1:50));
    title(['OdpowiedŸ dla systemu nr: ' num2str(i)]);
    xlabel('Czas [s]'); grid;
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
