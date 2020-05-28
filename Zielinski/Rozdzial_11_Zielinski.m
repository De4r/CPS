% Rodzial 11 Zielniski
%                       Mateusz Krupnik

% Filtry cyfrowe na podstawie filtrów analogowych z rozdzia³u 6.
clc; close all; clear all;
% Analogowe filtry Butterwortha i Czebyszewa
apass = 3;  % nielinowoœæ pasma przepustowego w dB
astop = 60; % t³umienie w paœmie zaprowym
fpr = 2000; % czestotliwoœæ probkowania
fmx = 1000; % maks. sk³adowa czêstotliwoœci fpr/2

%% Filtry Butterwortha Hp, BP i BS - filtry oparte na okreglu w lewej polpl.
typ = "Butterworth ";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = "LowPass ";
        fpass = 200; fstop = 300; % Hz
        fpass = fc2fa(fpass, fpr);
        fstop = fc2fa(fstop, fpr);
        ws = fstop/fpass;   % czestotliwoœæ znorm. s=s'/w0, w0=2*pi*fpass
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = "HighPass ";
        fstop = 700; % czês. pasma przepustowego odpowiadaj¹ca astop
        fpass = 800; % czês. pasma zaporowego odpowiadaj¹ca apass
        fpass = fc2fa(fpass, fpr);
        fstop = fc2fa(fstop, fpr);
        ws = fpass/fstop; % transformacja czês.i: s=w0/s', w0=2*pi*fpass
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = "BandPass ";
        fs1 = 300; % dolna czêstotliwoœæ stop
        fp1 = 400; % dolna czêstotliwoœæ pass
        fp2 = 600; % górna czêstotliwoœæ pass
        fs2 = 700; % górna czêstotliwoœæ stop
        % transformacja czêstotliwoœci
        fs1 = fc2fa(fs1, fpr);
        fp1 = fc2fa(fp1, fpr);
        fs2 = fc2fa(fs2, fpr);
        fp2 = fc2fa(fp2, fpr);
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    % filtr srodkowozaporowy BS
    else
        filtr = "BandStop ";
        fp1 = 200; % dolna czêstotliwoœæ filtra pasmowego
        fs1 = 300; % dolna czêstotliwoœæ filtra pasmowego
        fs2 = 700; % górna czêstotliwoœæ filtra pasmowego
        fp2 = 800; % górna czêstotliwoœæ filtra pasmowego
        
        fs1 = fc2fa(fs1, fpr);
        fp1 = fc2fa(fp1, fpr);
        fs2 = fc2fa(fs2, fpr);
        fp2 = fc2fa(fp2, fpr);
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
    
    disp([num2str(i) ' - FILTR CYFROWY: ' filtr ' - ' typ]); 
    
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
    figure(7*i-6)
    plot( real(p), imag(p), 'x' ); grid;
    title(["Po³o¿enie biegunów dla filtra analogowego:" + filtr + typ]);
    xlabel('real'); ylabel('imag');
    
    % Porównanie z Matlabem
    [NN,ww0] = buttord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystykê czêstotliwoœciow¹ H(w)=B(w)/A(w)
    w = 0 : 0.005 : 2;	% zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
    H = freqs(b,a,w);	% alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(7*i-5); subplot(211);
    plot(w,abs(H)); grid; title(["Modu³ prototypu anaglowego LowPass dla" + filtr + typ]);
    xlabel('Pulsacja [rad/sek]');
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(["Modu³ prototypu analogowego LowPass w dB dla " + filtr + typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    
    % Transformata czêstotliwoœci filtra analogowego: prototyp unormowany --> wynikowy filtr
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end % LowPass to LowPass: s=s/w0
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end % LowPass to HighPass: s=w0/s
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    b=wzm*poly(z); a=poly(p);
    
    % Poka¿ zera i bieguny po transformacji czêstoliwoœci
    figure(7*i-4)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(["Po³o¿enie biegunów dla filtra analogowego" + filtr + typ]);
    xlabel('real'); ylabel('imag');

    printsys(b,a,'s');
    % Koñcowa charakterystyka czêstoliwoœciowa
    NF = 1000; % ile punktów
    fmin = 0; % dolna czêstotliwoœæ
    fmax = 5000; % górna czêstotliwoœæ
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie czêstotliwoœci
    w = 2*pi*f; % wszystkie pulasacje
    H = freqs(b,a,w); % alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(7*i-3); subplot(211);
    plot( f,abs(H), f_ps, wzm_ps,'ro');
    grid; title(["Modu³ dla filtra: " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(["Modu³ w dB filtra " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
    
    % Filtr cyfrowy
    [ zc, pc, wzmc ] = bilinearTZ(z, p, wzm, fpr);
    bc = wzmc*poly(zc); ac = poly(pc);
    
    % Wykres zer filtra cyfrowego
    NP = 1000; fi=2*pi*(0:1:NP-1)/NP; x=sin(fi); y=cos(fi);
    figure(7*i-2)
    plot(x,y,'-k',real(zc),imag(zc),'or',real(pc),imag(pc),'xb');
    title(["ZERA i BIEGUNY filtra cyfrowego " + filtr + typ]);
    grid;
    
    % Charakterystyka filtra cyfrowego 
    NF = 1000; fmin = 0; fmax = fmx; f = fmin : (fmax-fmin)/(NF-1) : fmax;
    w = 2*pi*f/fpr; H = freqz(bc,ac,w);
    Habs=abs(H); HdB=20*log10(Habs); Hfa=unwrap(angle(H));
    f_ps = (fpr/pi)*atan(pi*f_ps/fpr);
    
    figure(7*i-1); subplot(211);
    plot( f,Habs, f_ps, wzm_ps,'ro');
    grid; title(["Modu³ dla filtra cyfrowego: " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f, HdB, f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(["Modu³ w dB filtra cyfrowego " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f, Hfa); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
    
    % Odpowiedz filtra na delte Kroneckera
    Nx=200; x = zeros(1,Nx); x(1)=1;
    M=length(bc); N=length(ac);
    ac=ac(2:N); N=N-1;
    bx=zeros(1,M); by=zeros(1,N); y=[];
    for n=1:Nx
        bx = [ x(n) bx(1:M-1)];
        y(n) = sum(bx .* bc) - sum(by .* ac);
        by = [ y(n) by(1:N-1) ];
    end
    n=0:Nx-1;
    figure(7*i);
    plot(n,y); grid;
    title(["Odp. impulsowa filtra cyfrowego " + filtr + typ]);
    xlabel('Probki - n'); ylabel('Ampl.');
end

%% Filtry na podstawie prototypu Czebyszewa I typu
typ = "Czebyszewa I ";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = "LowPass ";
        fpass = 200; fstop = 300; % Hz
        fpass = fc2fa(fpass, fpr);
        fstop = fc2fa(fstop, fpr);
        ws = fstop/fpass;   % czestotliwoœæ znorm. s=s'/w0, w0=2*pi*fpass
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = "HighPass ";
        fstop = 700; % czês. pasma przepustowego odpowiadaj¹ca astop
        fpass = 800; % czês. pasma zaporowego odpowiadaj¹ca apass
        fpass = fc2fa(fpass, fpr);
        fstop = fc2fa(fstop, fpr);
        ws = fpass/fstop; % transformacja czês.i: s=w0/s', w0=2*pi*fpass
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = "BandPass ";
        fs1 = 300; % dolna czêstotliwoœæ stop
        fp1 = 400; % dolna czêstotliwoœæ pass
        fp2 = 600; % górna czêstotliwoœæ pass
        fs2 = 700; % górna czêstotliwoœæ stop
        % transformacja czêstotliwoœci
        fs1 = fc2fa(fs1, fpr);
        fp1 = fc2fa(fp1, fpr);
        fs2 = fc2fa(fs2, fpr);
        fp2 = fc2fa(fp2, fpr);
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    % filtr srodkowozaporowy BS
    else
        filtr = "BandStop ";
        fp1 = 200; % dolna czêstotliwoœæ filtra pasmowego
        fs1 = 300; % dolna czêstotliwoœæ filtra pasmowego
        fs2 = 700; % górna czêstotliwoœæ filtra pasmowego
        fp2 = 800; % górna czêstotliwoœæ filtra pasmowego
        
        fs1 = fc2fa(fs1, fpr);
        fp1 = fc2fa(fp1, fpr);
        fs2 = fc2fa(fs2, fpr);
        fp2 = fc2fa(fp2, fpr);
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
    
    disp([num2str(i) ' - FILTR CYFROWY: ' filtr ' - ' typ]); 
    
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
    
    figure(7*i-6+28)
    plot( real(p), imag(p), 'x' ); grid;
    title(["Po³o¿enie biegunów dla filtra analogowego:" + filtr + typ]);
    xlabel('real'); ylabel('imag');
    
    % Porównanie z Matlabem
    [NN,ww0] = cheb1ord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystykê czêstotliwoœciow¹ H(w)=B(w)/A(w)
    w = 0 : 0.005 : 2;	% zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
    H = freqs(b,a,w);	% alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(7*i-5+28); subplot(211);
    plot(w,abs(H)); grid; title(["Modu³ prototypu anaglowego LowPass dla" + filtr + typ]);
    xlabel('Pulsacja [rad/sek]');
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(["Modu³ prototypu analogowego LowPass w dB dla " + filtr + typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    
    % Transformata czêstotliwoœci filtra analogowego: prototyp unormowany --> wynikowy filtr
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end % LowPass to LowPass: s=s/w0
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end % LowPass to HighPass: s=w0/s
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    b=wzm*poly(z); a=poly(p);
    
    % Poka¿ zera i bieguny po transformacji czêstoliwoœci
    figure(7*i-4+28)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(["Po³o¿enie biegunów dla filtra analogowego" + filtr + typ]);
    xlabel('real'); ylabel('imag');

    printsys(b,a,'s');
    % Koñcowa charakterystyka czêstoliwoœciowa
    NF = 1000; % ile punktów
    fmin = 0; % dolna czêstotliwoœæ
    fmax = 5000; % górna czêstotliwoœæ
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie czêstotliwoœci
    w = 2*pi*f; % wszystkie pulasacje
    H = freqs(b,a,w); % alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(7*i-3+28); subplot(211);
    plot( f,abs(H), f_ps, wzm_ps,'ro');
    grid; title(["Modu³ dla filtra: " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(["Modu³ w dB filtra " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
    
    % Filtr cyfrowy
    [ zc, pc, wzmc ] = bilinearTZ(z, p, wzm, fpr);
    bc = wzmc*poly(zc); ac = poly(pc);
    
    % Wykres zer filtra cyfrowego
    NP = 1000; fi=2*pi*(0:1:NP-1)/NP; x=sin(fi); y=cos(fi);
    figure(7*i-2+28)
    plot(x,y,'-k',real(zc),imag(zc),'or',real(pc),imag(pc),'xb');
    title(["ZERA i BIEGUNY filtra cyfrowego " + filtr + typ]);
    grid;
    
    % Charakterystyka filtra cyfrowego 
    NF = 1000; fmin = 0; fmax = fmx; f = fmin : (fmax-fmin)/(NF-1) : fmax;
    w = 2*pi*f/fpr; H = freqz(bc,ac,w);
    Habs=abs(H); HdB=20*log10(Habs); Hfa=unwrap(angle(H));
    f_ps = (fpr/pi)*atan(pi*f_ps/fpr);
    
    figure(7*i-1+28); subplot(211);
    plot( f,Habs, f_ps, wzm_ps,'ro');
    grid; title(["Modu³ dla filtra cyfrowego: " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f, HdB, f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(["Modu³ w dB filtra cyfrowego " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f, Hfa); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
    
    % Odpowiedz filtra na delte Kroneckera
    Nx=200; x = zeros(1,Nx); x(1)=1;
    M=length(bc); N=length(ac);
    ac=ac(2:N); N=N-1;
    bx=zeros(1,M); by=zeros(1,N); y=[];
    for n=1:Nx
        bx = [ x(n) bx(1:M-1)];
        y(n) = sum(bx .* bc) - sum(by .* ac);
        by = [ y(n) by(1:N-1) ];
    end
    n=0:Nx-1;
    figure(7*i+28);
    plot(n,y); grid;
    title(["Odp. impulsowa filtra cyfrowego " + filtr + typ]);
    xlabel('Probki - n'); ylabel('Ampl.');
end

%% Filtry na podstawie prototypu Czebyszewa II typu
typ = "Czebyszewa II ";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = "LowPass ";
        fpass = 200; fstop = 300; % Hz
        fpass = fc2fa(fpass, fpr);
        fstop = fc2fa(fstop, fpr);
        ws = fstop/fpass;   % czestotliwoœæ znorm. s=s'/w0, w0=2*pi*fpass
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = "HighPass ";
        fstop = 700; % czês. pasma przepustowego odpowiadaj¹ca astop
        fpass = 800; % czês. pasma zaporowego odpowiadaj¹ca apass
        fpass = fc2fa(fpass, fpr);
        fstop = fc2fa(fstop, fpr);
        ws = fpass/fstop; % transformacja czês.i: s=w0/s', w0=2*pi*fpass
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = "BandPass ";
        fs1 = 300; % dolna czêstotliwoœæ stop
        fp1 = 400; % dolna czêstotliwoœæ pass
        fp2 = 600; % górna czêstotliwoœæ pass
        fs2 = 700; % górna czêstotliwoœæ stop
        % transformacja czêstotliwoœci
        fs1 = fc2fa(fs1, fpr);
        fp1 = fc2fa(fp1, fpr);
        fs2 = fc2fa(fs2, fpr);
        fp2 = fc2fa(fp2, fpr);
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    % filtr srodkowozaporowy BS
    else
        filtr = "BandStop ";
        fp1 = 200; % dolna czêstotliwoœæ filtra pasmowego
        fs1 = 300; % dolna czêstotliwoœæ filtra pasmowego
        fs2 = 700; % górna czêstotliwoœæ filtra pasmowego
        fp2 = 800; % górna czêstotliwoœæ filtra pasmowego
        
        fs1 = fc2fa(fs1, fpr);
        fp1 = fc2fa(fp1, fpr);
        fs2 = fc2fa(fs2, fpr);
        fp2 = fc2fa(fp2, fpr);
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
    
    disp([num2str(i) ' - FILTR CYFROWY: ' filtr ' - ' typ]); 
    
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

    
    figure(7*i-6+56)
    plot( real(p), imag(p), 'x' ); grid;
    title(["Po³o¿enie biegunów dla filtra analogowego:" + filtr + typ]);
    xlabel('real'); ylabel('imag');
    
    % Porównanie z Matlabem
    [NN,ww0] = cheb2ord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystykê czêstotliwoœciow¹ H(w)=B(w)/A(w)
    w = 0 : 0.005 : 2;	% zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
    H = freqs(b,a,w);	% alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(7*i-5+56); subplot(211);
    plot(w,abs(H)); grid; title(["Modu³ prototypu anaglowego LowPass dla" + filtr + typ]);
    xlabel('Pulsacja [rad/sek]');
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(["Modu³ prototypu analogowego LowPass w dB dla " + filtr + typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    
    % Transformata czêstotliwoœci filtra analogowego: prototyp unormowany --> wynikowy filtr
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end % LowPass to LowPass: s=s/w0
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end % LowPass to HighPass: s=w0/s
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    b=wzm*poly(z); a=poly(p);
    
    % Poka¿ zera i bieguny po transformacji czêstoliwoœci
    figure(7*i-4+56)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(["Po³o¿enie biegunów dla filtra analogowego" + filtr + typ]);
    xlabel('real'); ylabel('imag');

    printsys(b,a,'s');
    % Koñcowa charakterystyka czêstoliwoœciowa
    NF = 1000; % ile punktów
    fmin = 0; % dolna czêstotliwoœæ
    fmax = 5000; % górna czêstotliwoœæ
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie czêstotliwoœci
    w = 2*pi*f; % wszystkie pulasacje
    H = freqs(b,a,w); % alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(7*i-3+56); subplot(211);
    plot( f,abs(H), f_ps, wzm_ps,'ro');
    grid; title(["Modu³ dla filtra: " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(["Modu³ w dB filtra " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
    
    % Filtr cyfrowy
    [ zc, pc, wzmc ] = bilinearTZ(z, p, wzm, fpr);
    bc = wzmc*poly(zc); ac = poly(pc);
    
    % Wykres zer filtra cyfrowego
    NP = 1000; fi=2*pi*(0:1:NP-1)/NP; x=sin(fi); y=cos(fi);
    figure(7*i-2+56)
    plot(x,y,'-k',real(zc),imag(zc),'or',real(pc),imag(pc),'xb');
    title(["ZERA i BIEGUNY filtra cyfrowego " + filtr + typ]);
    grid;
    
    % Charakterystyka filtra cyfrowego 
    NF = 1000; fmin = 0; fmax = fmx; f = fmin : (fmax-fmin)/(NF-1) : fmax;
    w = 2*pi*f/fpr; H = freqz(bc,ac,w);
    Habs=abs(H); HdB=20*log10(Habs); Hfa=unwrap(angle(H));
    f_ps = (fpr/pi)*atan(pi*f_ps/fpr);
    
    figure(7*i-1+56); subplot(211);
    plot( f,Habs, f_ps, wzm_ps,'ro');
    grid; title(["Modu³ dla filtra cyfrowego: " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]');
    subplot(212);
    plot(f, HdB, f_ps, wzmdB_ps,'ro'); axis([fmin,fmax,-100,20]);
    grid; title(["Modu³ w dB filtra cyfrowego " + filtr + typ]);
    xlabel('Czestotliwoœæ [Hz]'); ylabel('dB');
    plot(f, Hfa); grid; title('Faza');
    xlabel('Czestotliwoœæ [Hz]'); ylabel('[rad]');
    
    % Odpowiedz filtra na delte Kroneckera
    Nx=200; x = zeros(1,Nx); x(1)=1;
    M=length(bc); N=length(ac);
    ac=ac(2:N); N=N-1;
    bx=zeros(1,M); by=zeros(1,N); y=[];
    for n=1:Nx
        bx = [ x(n) bx(1:M-1)];
        y(n) = sum(bx .* bc) - sum(by .* ac);
        by = [ y(n) by(1:N-1) ];
    end
    n=0:Nx-1;
    figure(7*i+56);
    plot(n,y); grid;
    title(["Odp. impulsowa filtra cyfrowego " + filtr + typ]);
    xlabel('Probki - n'); ylabel('Ampl.');
end

%%%%% DEFINICJE FUNKCJI %%%%%%%%

% Transformacja biliniowa filtru analogowego do cyfrowego
function [zz, pp, wzm] = bilinearTZ(z, p , wzm, fpr)
    pp = []; zz = [];
    for k=1:length(z)
       zz = [ zz (2*fpr+z(k))/(2*fpr-z(k)) ];
       wzm = wzm*(2*fpr-z(k));
    end
    for k=1:length(p)
        pp = [ pp (2*fpr+p(k))/(2*fpr-p(k)) ];
        wzm = wzm/(2*fpr-p(k));
    end
    l1 = length(p) - length(z);
    l2 = length(z) - length(p);
    if (l1 > 0) zz = [ zz -1*ones(1, l1) ]; end
    if (l2 > 0) pp = [ pp -1*ones(1, l2) ]; end
end

% PRzeliczenie czestotliwosci cyfrowej na analogowa
function fn = fc2fa(f, fpr)
    fn = 2*fpr*tan(pi*f/fpr)/(2*pi);
end

%%% Funkcje z rozdz 6.

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