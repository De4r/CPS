% Rodzial 5 Zielniski
%                       Mateusz Krupnik

% Projektowanie filtrów metoda zer i biegunow
clc; clear all; close all;
% Przyk³ad 1: projekt filtra pasmowoprzepustowego
% o wpass1 = 9.5 rd, wpass2 = 10.5 rd

z1 = 5; z2 = 15;                % ZERA na osi urojonej
z = j*[ -z2, -z1, z1, z2 ];
odl = 0.5; p1 = 9.5; p2 = 10.5;	% BIEGUNY w pobli¿u osi urojonej
p = [ -odl-j*p2, -odl-j*p1, -odl+j*p1, -odl+j*p2 ];
WMAX=20; TMAX=20;               % max pulsacja, max czas obserwacji
show_results(z, p, WMAX, TMAX, "przyklad 1.");

%% Przyk³ad 2: znajdowanie zer i biegunów zadanej transmitancji
b=[ 0.66667 0 1 ];              % wspó³czynniki licznika transmitancji
a=[ 4.0001 5.0081 3.1650 1 ];   % wspó³czynniki mianownika transmitancji
[z,p,wzm] = tf2zp(b,a);         % wspó³czynniki wielomianów -> zera wielomianów
z = z'; p = p';                 % wektor pionowy -> wektor poziomy
WMAX=5; TMAX=25;                % max pulsacja, max czas obserwacji
show_results(z, p, WMAX, TMAX, "przyklad 2.");

%% Przyk³ad 3: projekt filtra górnoprzepustowego
z1 = 0;                 % ZERA na osi urojonej
z2 = 0+j*1; z3 = 0-j*1;
z4 = 0+j*2; z5 = 0-j*2;
z6 = 0+j*3; z7 = 0-j*3;
z = [ z1 z2 z3 z4 z5 z6 z7 ];
p1 = -1;                % BIEGUNY w pobli¿u osi urojonej
p2 = -1+j*1; p3 = -1-j*1;
p4 = -1+j*2; p5 = -1-j*2;
p6 = -1+j*3; p7 = -1-j*3;
p = [ p1 p2 p3 p4 p5 p6 p7 ];
WMAX=20; TMAX=5;        % max pulsacja, max czas obserwacji
show_results(z, p, WMAX, TMAX, "przyklad 3.");



%%%%%%%%%%%%%%%% FUNKCJA PREZENTACJI %%%%%%%%%%%%%%%%%%%%%%%%%
function show_results(z, p, WMAX, TMAX, title_text)
    figure();
    plot(real(z), imag(z), 'or', real(p), imag(p), 'xb'); grid;
    title(['Zera (o) i bieguny (x):' title_text]);
    xlabel('Real'); ylabel('Imag [rd/s]');
    w = 0 : 0.01 : WMAX; % wybrane pulsacje widma
    [b,a] = zp2tf(z',p',1); % zera, bieguny -> wspólczynniki wielomianów
    H = freqs(b,a,w); % wyznaczenie widma transmitancji dla zadanego w
    Hm = abs(H); HmdB = 20*log10(Hm); % modu³ transmitancji
    Hf = angle(H); Hfu = unwrap(Hf); % faza transmitancji
    figure(); sgtitle(['Charakterystyki fazowe dla: ' title_text]);
    subplot(221);
    plot(w,Hm,'k'); grid;
    title('Ch-ka amplitudowa'); xlabel('Czestoœæ w[rd/s]');
    subplot(222);
    plot(w,HmdB,'k'); grid;
    title('Ch-ka amplitudowa w dB'); xlabel('Czestoœæ w[rd/s]');
    subplot(223);
    plot(w,Hf,'k'); grid;
    title('Ch-ka fazowa'); xlabel('Czestoœæ w[rd/s]'); ylabel('[rd]');
    subplot(224);
    plot(w,Hfu,'k'); grid;
    title('Ch-ka fazowa unwrap'); xlabel('Czestoœæ w[rd/s]'); ylabel('[rd]');

    % OdpowiedŸ impulsowa
    h = impulse(b,a,TMAX); % funkcja z przybornika CONTROL
    dt = TMAX/(length(h)-1); th = 0 : dt : TMAX;
    figure(); sgtitle(['Charakterystyki czasowe dla: ' title_text]);
    subplot(211);
    plot(th,h,'k'); grid; title('OdpowiedŸ impulsowa'); xlabel('Czas t[s]');
    % OdpowiedŸ na skok jednostkowy
    u = step(b,a,TMAX); % funkcja z przybornika CONTROL
    dt = TMAX/(length(u)-1); tu = 0 : dt : TMAX;
    subplot(212);
    plot(tu,u,'k'); grid; title('OdpowiedŸ skokowa'); xlabel('Czas t[s]');
end