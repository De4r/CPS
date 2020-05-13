% Rodzial 6 Zielniski
%                       Mateusz Krupnik

clc; close all; clear all;
% Analogowe filtry Butterwortha i Czebyszewa
apass = 1;  % nielinowoœæ pasma przepustowego w dB
astop = 50; % t³umienie w paœmie zaprowym

%% Filtry Butterwortha Hp, BP i BS - filtry oparte na okreglu w lewej polpl.
typ = "Butterworth";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = 'LowPass';
        fpass = 1000; fstop = 4000; % Hz
        ws = fstop/fpass;   % czestotliwoœæ znorm. s=s'/w0, w0=2*pi*fpass
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = 'HighPass';
        fstop = 2000; % czês. pasma przepustowego odpowiadaj¹ca astop
        fpass = 3000; % czês. pasma zaporowego odpowiadaj¹ca apass
        ws = fpass/fstop; % transformacja czês.i: s=w0/s', w0=2*pi*fpass
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = 'BandPass';
        fs1 = 1500; % dolna czêstotliwoœæ stop
        fp1 = 2000; % dolna czêstotliwoœæ pass
        fp2 = 3000; % górna czêstotliwoœæ pass
        fs2 = 3500; % górna czêstotliwoœæ stop
        % transformacja czêstotliwoœci
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    % filtr srodkowozaporowy BS
    else
        filtr = 'BandStop';
        fp1 = 1500; % dolna czêstotliwoœæ filtra pasmowego
        fs1 = 2000; % dolna czêstotliwoœæ filtra pasmowego
        fs2 = 3000; % górna czêstotliwoœæ filtra pasmowego
        fp2 = 3500; % górna czêstotliwoœæ filtra pasmowego
        % transformacja czêstotliwoœci
        ws1t = (fs1*(fp2-fp1)) / (fs1^2 - fp1*fp2); 
        ws2t = (fs2*(fp2-fp1)) / (fs2^2 - fp1*fp2);
        ws = min(abs(ws1t), abs(ws2t));
    end
    % Przeliczenie na wartoœæ bezwzglêdn¹
    wzm_p = 10^(-apass/20);
    wzm_s = 10^(-astop/20);
    
    if ( (i==1) || (i==2))
        vp = 2*pi*fpass;
        vs = 2*pi*fstop;
        f_ps = [fpass, fstop];
        wzm_ps = [wzm_p, wzm_s];
        wzmdB_ps = [-apass, -astop];
    else
        vp = 2*pi*[ fp1 fp2 ];
        vs = 2*pi*[ fs1 fs2 ];
        vc = 2*pi*sqrt(fp1*fp2);	% pulsacja œrodka
        % szerokoœæ filtra wokó³ vc
        f_ps = [fp1,fp2,fs1,fs2];
        dv = 2*pi*(fp2-fp1);
        wzm_ps = [wzm_p, wzm_p, wzm_s, wzm_s];
        wzmdB_ps = [-apass, -apass, -astop, -astop];
    end
    
    disp([num2str(i) ' - FILTR: ' filtr ' - ' typ]); 
    
    % Obliczenie parametrów N i w0 - funckja buttord
    wp = 1;
    N = ceil(log10( (10^(astop/10)-1) /(10^(apass/10)-1) ) / (2*log10(ws/wp)) );
    w0 = ws / (10^(astop/10)-1)^(1/(2*N));

    % Obliczenie biegunów trans. prototypu - funcja buttap i zp2tf
    dfi0 = (2*pi)/(2*N);                    % k¹t „kawa³ka tortu”
    fi = pi/2 + dfi0/2 + (0 : N-1)*dfi0;    % k¹ty biegunów
    p = w0*exp(j*fi);                       % bieguny
    z = [];                                 % zera
    wzm = real( prod(-p) );                 % wzmocnienie
    a = poly(p);                % bieguny --> wsp wielomianu mianownika A(z)
    b = wzm;                    % wielomian licznika B(z)
%     z, p, b, a
    figure(4*i-3)
    plot( real(p), imag(p), 'x' ); grid;
    title(['Po³o¿enie biegunów dla filtra:' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
    
    % Porównanie z Matlabem
    [NN,ww0] = buttord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystykê czêstotliwoœciow¹ H(w)=B(w)/A(w)
    w = 0 : 0.005 : 2;	% zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
    H = freqs(b,a,w);	% alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i-2); subplot(211);
    plot(w,abs(H)); grid; title(['Modu³ prototypu LowPass dla' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]');
    
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(['Modu³ prototypu LowPass w dB dla ' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    % Transformata czêstotliwoœci filtra analogowego: prototyp unormowany --> wynikowy filtr
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end % LowPass to LowPass: s=s/w0
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end % LowPass to HighPass: s=w0/s
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    b=wzm*poly(z); a=poly(p);
    % Poka¿ zera i bieguny po transformacji czêstoliwoœci
    figure(4*i-1)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(['Po³o¿enie biegunów dla filtra ' filtr ' - ' typ]); xlabel('real'); ylabel('imag');
%     p, z
%     a, b
    printsys(b,a,'s');
    % Koñcowa charakterystyka czêstoliwoœciowa
    NF = 1000; % ile punktów
    fmin = 0; % dolna czêstotliwoœæ
    fmax = 5000; % górna czêstotliwoœæ
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie czêstotliwoœci
    w = 2*pi*f; % wszystkie pulasacje
    H = freqs(b,a,w); % alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i); subplot(211);
    plot( f,abs(H), f_ps, wzm_ps,'ro');
    grid; title(['Modu³ dla filtra: ' filtr ' - ' typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(['Modu³ w dB filtra ' filtr ' - ' typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
end

%% Filtry na podstawie prototypu Czebyszewa I typu
typ = "Czebyszewa I";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = 'LowPass';
        fpass = 1000; fstop = 4000; % Hz
        ws = fstop/fpass;   % czestotliwoœæ znorm. s=s'/w0, w0=2*pi*fpass
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = 'HighPass';
        fstop = 2000; % czês. pasma przepustowego odpowiadaj¹ca astop
        fpass = 3000; % czês. pasma zaporowego odpowiadaj¹ca apass
        ws = fpass/fstop; % transformacja czês.i: s=w0/s', w0=2*pi*fpass
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = 'BandPass';
        fs1 = 1500; % dolna czêstotliwoœæ stop
        fp1 = 2000; % dolna czêstotliwoœæ pass
        fp2 = 3000; % górna czêstotliwoœæ pass
        fs2 = 3500; % górna czêstotliwoœæ stop
        % transformacja czêstotliwoœci
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    % filtr srodkowozaporowy BS
    else
        filtr = 'BandStop';
        fp1 = 1500; % dolna czêstotliwoœæ filtra pasmowego
        fs1 = 2000; % dolna czêstotliwoœæ filtra pasmowego
        fs2 = 3000; % górna czêstotliwoœæ filtra pasmowego
        fp2 = 3500; % górna czêstotliwoœæ filtra pasmowego
        % transformacja czêstotliwoœci
        ws1t = (fs1*(fp2-fp1)) / (fs1^2 - fp1*fp2); 
        ws2t = (fs2*(fp2-fp1)) / (fs2^2 - fp1*fp2);
        ws = min(abs(ws1t), abs(ws2t));
    end
    % Przeliczenie na wartoœæ bezwzglêdn¹
    wzm_p = 10^(-apass/20);
    wzm_s = 10^(-astop/20);
    
    if ( (i==1) || (i==2))
        vp = 2*pi*fpass;
        vs = 2*pi*fstop;
        f_ps = [fpass, fstop];
        wzm_ps = [wzm_p, wzm_s];
        wzmdB_ps = [-apass, -astop];
    else
        vp = 2*pi*[ fp1 fp2 ];
        vs = 2*pi*[ fs1 fs2 ];
        vc = 2*pi*sqrt(fp1*fp2);	% pulsacja œrodka
        % szerokoœæ filtra wokó³ vc
        f_ps = [fp1,fp2,fs1,fs2];
        dv = 2*pi*(fp2-fp1);
        wzm_ps = [wzm_p, wzm_p, wzm_s, wzm_s];
        wzmdB_ps = [-apass, -apass, -astop, -astop];
    end
    
    disp([num2str(i) ' - FILTR: ' filtr ' - ' typ]); 
    
    % Obliczenie parametrów N i w0 
    wp = 1;
    Nreal = acosh(sqrt((10^(astop/10)-1) / (10^(apass/10)-1))) / acosh(ws/wp);
    N = ceil(Nreal);
    epsi = sqrt(10^(apass/10)-1);
    D = asinh(1/epsi)/N; R1 = sinh(D); R2 = cosh(D);

    % Obliczenie biegunów trans. prototypu - funcja buttap i zp2tf
    dfi0 = (2*pi)/(2*N);                    % k¹t „kawa³ka tortu”
    fi = pi/2 + dfi0/2 + (0 : N-1)*dfi0;    % k¹ty biegunów
    p1 = R1 * exp(j*fi);                    % bieguny R1
    p2 = R2 * exp(j*fi);                    % bieguny R2
    p = real(p1) +1i*imag(p2);               % Polaczone bieguny
    z = [];                                 % zera
    wzm = prod(-p);                 % wzmocnienie
    a = poly(p);                % bieguny --> wsp wielomianu mianownika A(z)
    b = wzm;                    % wielomian licznika B(z)
    if (rem(N,2)==0) b = b*10^(-apass/20); end

    figure(4*i-3+16)
    plot( real(p), imag(p), 'x' ); grid;
    title(['Po³o¿enie biegunów dla filtra:' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
    
    % Porównanie z Matlabem
    [NN,ww0] = cheb1ord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystykê czêstotliwoœciow¹ H(w)=B(w)/A(w)
    w = 0 : 0.005 : 2;	% zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
    H = freqs(b,a,w);	% alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i-2+16); subplot(211);
    plot(w,abs(H)); grid; title(['Modu³ prototypu LowPass dla' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]');
    
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(['Modu³ prototypu LowPass w dB dla ' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    % Transformata czêstotliwoœci filtra analogowego: prototyp unormowany --> wynikowy filtr
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end % LowPass to LowPass: s=s/w0
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end % LowPass to HighPass: s=w0/s
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    b=wzm*poly(z); a=poly(p);
    % Poka¿ zera i bieguny po transformacji czêstoliwoœci
    figure(4*i-1+16)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(['Po³o¿enie biegunów dla filtra ' filtr ' - ' typ]); xlabel('real'); ylabel('imag');
%     p, z
%     a, b
    printsys(b,a,'s');
    % Koñcowa charakterystyka czêstoliwoœciowa
    NF = 1000; % ile punktów
    fmin = 0; % dolna czêstotliwoœæ
    fmax = 5000; % górna czêstotliwoœæ
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie czêstotliwoœci
    w = 2*pi*f; % wszystkie pulasacje
    H = freqs(b,a,w); % alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i+16); subplot(211);
    plot(f, abs(H), f_ps, wzm_ps,'ro');
    grid; title(['Modu³ dla filtra: ' filtr ' - ' typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(['Modu³ w dB filtra ' filtr ' - ' typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
end

%% Filtry na podstawie prototypu Czebyszewa II typu
typ = "Czebyszewa II";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = 'LowPass';
        fpass = 1000; fstop = 4000; % Hz
        ws = fstop/fpass;   % czestotliwoœæ znorm. s=s'/w0, w0=2*pi*fpass
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = 'HighPass';
        fstop = 2000; % czês. pasma przepustowego odpowiadaj¹ca astop
        fpass = 3000; % czês. pasma zaporowego odpowiadaj¹ca apass
        ws = fpass/fstop; % transformacja czês.i: s=w0/s', w0=2*pi*fpass
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = 'BandPass';
        fs1 = 1500; % dolna czêstotliwoœæ stop
        fp1 = 2000; % dolna czêstotliwoœæ pass
        fp2 = 3000; % górna czêstotliwoœæ pass
        fs2 = 3500; % górna czêstotliwoœæ stop
        % transformacja czêstotliwoœci
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    % filtr srodkowozaporowy BS
    else
        filtr = 'BandStop';
        fp1 = 1500; % dolna czêstotliwoœæ filtra pasmowego
        fs1 = 2000; % dolna czêstotliwoœæ filtra pasmowego
        fs2 = 3000; % górna czêstotliwoœæ filtra pasmowego
        fp2 = 3500; % górna czêstotliwoœæ filtra pasmowego
        % transformacja czêstotliwoœci
        ws1t = (fs1*(fp2-fp1)) / (fs1^2 - fp1*fp2); 
        ws2t = (fs2*(fp2-fp1)) / (fs2^2 - fp1*fp2);
        ws = min(abs(ws1t), abs(ws2t));
    end
    % Przeliczenie na wartoœæ bezwzglêdn¹
    wzm_p = 10^(-apass/20);
    wzm_s = 10^(-astop/20);
    
    if ( (i==1) || (i==2))
        vp = 2*pi*fpass;
        vs = 2*pi*fstop;
        f_ps = [fpass, fstop];
        wzm_ps = [wzm_p, wzm_s];
        wzmdB_ps = [-apass, -astop];
    else
        vp = 2*pi*[ fp1 fp2 ];
        vs = 2*pi*[ fs1 fs2 ];
        vc = 2*pi*sqrt(fp1*fp2);	% pulsacja œrodka
        % szerokoœæ filtra wokó³ vc
        f_ps = [fp1,fp2,fs1,fs2];
        dv = 2*pi*(fp2-fp1);
        wzm_ps = [wzm_p, wzm_p, wzm_s, wzm_s];
        wzmdB_ps = [-apass, -apass, -astop, -astop];
    end
    
    disp([num2str(i) ' - FILTR: ' filtr ' - ' typ]); 
    
    % Obliczenie parametrów N i w0
    wp = 1;
    Nreal = acosh(sqrt((10^(astop/10)-1) / (10^(apass/10)-1))) / acosh(ws/wp);
    N = ceil(Nreal);
    epsi = sqrt(1 / (10^(apass/10)-1));
    D = asinh(1/epsi)/N; R1 = sinh(D); R2 = cosh(D);

    % Obliczenie biegunów trans. prototypu - funcja buttap i zp2tf
    dfi0 = (2*pi)/(2*N);                    % k¹t „kawa³ka tortu”
    fi = pi/2 + dfi0/2 + (0 : N-1)*dfi0;    % k¹ty biegunów
    p1 = R1 * exp(j*fi);                    % bieguny R1
    p2 = R2 * exp(j*fi);                    % bieguny R2
    p = real(p1) +1i*imag(p2);               % Polaczone bieguny
    z = 1i*sin(fi)                                 % zera
    wzm = prod(-z) / prod(-p);                 % wzmocnienie
    z = 1./z; p = 1./p;
    a = poly(p);                % bieguny --> wsp wielomianu mianownika A(z)
    b = wzm*poly(z);                    % wielomian licznika B(z)

    figure(4*i-3+32)
    plot( real(p), imag(p), 'x' ); grid;
    title(['Po³o¿enie biegunów dla filtra:' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
    
    % Porównanie z Matlabem
    [NN,ww0] = cheb2ord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystykê czêstotliwoœciow¹ H(w)=B(w)/A(w)
    w = 0 : 0.005 : 2;	% zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
    H = freqs(b,a,w);	% alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i-2+32); subplot(211);
    plot(w,abs(H)); grid; title(['Modu³ prototypu LowPass dla' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]');
    
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(['Modu³ prototypu LowPass w dB dla ' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    % Transformata czêstotliwoœci filtra analogowego: prototyp unormowany --> wynikowy filtr
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end % LowPass to LowPass: s=s/w0
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end % LowPass to HighPass: s=w0/s
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    b=wzm*poly(z); a=poly(p);
    % Poka¿ zera i bieguny po transformacji czêstoliwoœci
    figure(4*i-1+32)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(['Po³o¿enie biegunów dla filtra ' filtr ' - ' typ]); xlabel('real'); ylabel('imag');
%     p, z
%     a, b
    printsys(b,a,'s');
    % Koñcowa charakterystyka czêstoliwoœciowa
    NF = 1000; % ile punktów
    fmin = 0; % dolna czêstotliwoœæ
    fmax = 5000; % górna czêstotliwoœæ
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie czêstotliwoœci
    w = 2*pi*f; % wszystkie pulasacje
    H = freqs(b,a,w); % alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i+32); subplot(211);
    plot(f, abs(H), f_ps, wzm_ps,'ro');
    grid; title(['Modu³ dla filtra: ' filtr ' - ' typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(['Modu³ w dB filtra ' filtr ' - ' typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
end

%% Zaprojektowanie uk³adu elektronicznego dolnoprzepustowego filtra Butterwortha
% Dane projektowe
fpass = 8000;       % czêstotliwoœæ pasma przepustowego odpowiadaj¹ca apass
fstop = 22050;      % czêstotliwoœæ pasma zaporowego odpowiadaj¹ca astop
apass = 2;          % nieliniowoœæ pasma przepustowego w dB („zwis”)
astop = 40;         % t³umienie w paœmie zaporowym

wzm_p = 10^(-apass/20);	% t³umienie pass -> wzmocnienie pass
wzm_s = 10^(-astop/20);	% t³umienie stop -> wzmocnienie stop
ws = fstop/fpass;       % transformacja czêstotliwoœci: s=s'/w0, w0=2*pi*fpass
vp = 2*pi*fpass; vs = 2*pi*fstop;
f_ps = [fpass, fstop]; wzm_ps = [wzm_p, wzm_s]; wzmdB_ps = [-apass, -astop];

wp = 1;
Nreal = log10( (10^(astop/10)-1) / (10^(apass/10)-1) ) / (2*log10(ws/wp));
N = ceil(Nreal);
w0 = ws / (10^(astop/10)-1)^(1/(2*N));

dfi0 = (2*pi)/(2*N);                    % k¹t „kawa³ka tortu”
fi = pi/2 + dfi0/2 + (0 : N-1)*dfi0;	% k¹ty biegunów
p = w0*exp(1i*fi);                       % bieguny
z = [];                                 % zera
wzm = real(prod(-p));                   % wzmocnienie

figure(49)
plot( real(p), imag(p), 'x' ); grid;
title(['Po³o¿enie biegunów dla filtra dolnoprzepustowego']);
xlabel('real'); ylabel('imag');

b = wzm;                % wielomian licznika B(z)
a = poly(p);            % bieguny --> wsp wielomianu mianownika A(z)
printsys(b,a,'s');

% Porównanie z Matlabem
[NN,ww0] = buttord( vp, vs, apass, astop, 's' );
blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);

% Oblicz charakterystykê czêstotliwoœciow¹ H(w)=B(w)/A(w)
w = 0 : 0.005 : 2;	% zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
H = freqs(b,a,w);	% alternatywa: H = polyval( b,j*w)./polyval(a,j*w);

figure(50); subplot(211);
plot(w,abs(H)); grid; title(['Modu³ prototypu LowPass']);
xlabel('Pulsacja [rad/sek]');

subplot(212); plot(w,20*log10(abs(H))); grid;
title(['Modu³ prototypu LowPass w dB']);
xlabel('Pulsacja [rad/sek]'); ylabel('dB');

% Transformata czêstotliwoœci filtra analogowego: prototyp unormowany --> wynikowy filtr
[z,p,wzm] = lp2lpTZ(z,p,wzm,vp);     % LowPass to LowPass: s=s/w0

b=wzm*poly(z); a=poly(p);
% Poka¿ zera i bieguny po transformacji czêstoliwoœci
figure(51)
plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
title(['Po³o¿enie biegunów dla filtra LowPass']); xlabel('real'); ylabel('imag');
printsys(b,a,'s');

% Koñcowa charakterystyka czêstoliwoœciowa
NF = 1000; % ile punktów
fmin = 0; % dolna czêstotliwoœæ
fmax = 50000; % górna czêstotliwoœæ
f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie czêstotliwoœci
w = 2*pi*f; % wszystkie pulasacje
H = freqs(b,a,w); % alternatywa: H = polyval( b,j*w)./polyval(a,j*w);

figure(52); subplot(211);
plot( f,abs(H), f_ps, wzm_ps,'ro');
grid; title(['Modu³ dla filtra lowpass']);
xlabel('Czestotliwoœæ [Hz]');
subplot(212);
plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
grid; title(['Modu³ w dB filtra lowpass']);
xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
plot(f,unwrap(angle(H))); grid; title('Faza');
xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');

% Oblicz elementy uk³adu ze wzmacniaczami operacyjnymi
p1 = [ p(1) conj(p(1)) ];
p2 = [ p(2) conj(p(2)) ];
p3 = p(4);
aw1 = poly(p1), aw2 = poly(p2), aw3 = poly(p3);
C = 10^(-9); RA=10^4; Rwy = 10^4;

disp('=== Uk³ad 1===')
a = aw1;
a2=a(1); a1=a(2); a0=a(3);
R = 1/(C*sqrt(a0))
RB = (2-a1/sqrt(a0)) * RA
K1 = 1+RB/RA

disp('=== Uk³ad 2 ===')
a = aw2;
a2=a(1); a1=a(2); a0=a(3);
R = 1/(C*sqrt(a0))
RB = (2-a1/sqrt(a0)) * RA
K2 = 1+RB/RA

disp('=== Uk³ad 3 ===')
a = aw3;
a1=a(1); a0=a(2);
R=1/(C*a0)
K3=1

disp('=== Obci¹¿enie ===')
K=K1*K2*K3
G=1
Rx = (K/G)*Rwy
Ry = (G/K)/(1-G/K)*Rx



%%%%%%%%%%%%% DEFINICJE FUNKCJI ZAGNIEZDZONYCH %%%%%%%%%%%%%%

% Funkcja transformuj¹ca filtr LP znromalizowany na wymagany LP
function [zz,pp,wzm] = lp2lpTZ(z,p,wzm,w0)
    % LowPass to LowPass TZ
    zz = []; pp = [];
    for k=1:length(z)
        zz = [ zz z(k)*w0 ];
        wzm = wzm/w0;
    end
    for k=1:length(p)
        pp = [ pp p(k)*w0 ];
        wzm = wzm*w0;
    end
end

% Funkcja transformuj¹ca filtr LP znromalizowany na wymagany HP
function [zz,pp,wzm] = lp2hpTZ(z,p,wzm,w0)
    % LowPass to HighPass TZ
    zz = []; pp = [];
    for k=1:length(z)
        zz = [ zz w0/z(k) ];
        wzm = wzm*(-z(k));
    end
    for k=1:length(p)
        pp = [ pp w0/p(k) ];
        wzm = wzm/(-p(k));
    end
    for k=1:(length(p)-length(z))
        zz = [ zz 0 ];
    end
end


% Funkcja transformuj¹ca filtr LP znromalizowany na wymagany BP
function [zz,pp,wzm] = lp2bpTZ(z,p,wzm,w0,dw)
    % LowPass to BandPass TZ
    pp = []; zz = [];
    for k=1:length(z)
        zz = [ zz roots([ 1 -z(k)*dw w0^2])' ];
    wzm = wzm/dw;
    end
    for k=1:length(p)
        pp = [ pp roots([ 1 -p(k)*dw w0^2])' ];
        wzm = wzm*dw;
    end
    for k=1:(length(p)-length(z))
        zz = [ zz 0 ];
    end
end


% Funkcja transformuj¹ca filtr LP znromalizowany na wymagany BS
function [zz,pp,wzm] = lp2bsTZ(z,p,wzm,w0,dw)
    % LowPass to BandStop TZ
    zz = []; pp = [];
    for k=1:length(z)
        zz = [ zz roots([ 1 -dw/z(k) w0^2 ])' ];
        wzm = wzm*(-z(k));
    end
    for k=1:length(p)
        pp = [ pp roots([ 1 -dw/p(k) w0^2 ])' ];
        wzm = wzm/(-p(k));
    end
    for k=1:(length(p)-length(z))
        zz = [ zz roots([ 1 0 w0^2 ])' ];
    end
end


    