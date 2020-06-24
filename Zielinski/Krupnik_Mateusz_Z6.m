% Rodzial 6 Zielinski
%                       Mateusz Krupnik

clc; close all; clear all;
% Analogowe filtry Butterwortha i Czebyszewa
apass = 1;  % nielinowo�� pasma przepustowego w dB
astop = 50; % t�umienie w pa�mie zaprowym

%% Filtry Butterwortha HP, BP, BS - filtry oparte na okreglu w lewej polpl.
typ = "Butterworth";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = 'LowPass';
        fpass = 1000; fstop = 4000; % Hz
        ws = fstop/fpass;   % czestotliwo�� znorm. s=s'/w0, w0=2*pi*fpass
    
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = 'HighPass';
        fstop = 2000; % cz�s. pasma przepustowego odpowiadaj�ca astop
        fpass = 3000; % cz�s. pasma zaporowego odpowiadaj�ca apass
        ws = fpass/fstop; % transformacja cz�s.i: s=w0/s', w0=2*pi*fpass
    
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = 'BandPass';
        fs1 = 1500; % dolna cz�stotliwo�� stop
        fp1 = 2000; % dolna cz�stotliwo�� pass
        fp2 = 3000; % g�rna cz�stotliwo�� pass
        fs2 = 3500; % g�rna cz�stotliwo�� stop
        % transformacja cz�stotliwo�ci
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    
    % filtr srodkowozaporowy BS
    else
        filtr = 'BandStop';
        fp1 = 1500; % dolna cz�stotliwo�� filtra pasmowego
        fs1 = 2000; % dolna cz�stotliwo�� filtra pasmowego
        fs2 = 3000; % g�rna cz�stotliwo�� filtra pasmowego
        fp2 = 3500; % g�rna cz�stotliwo�� filtra pasmowego
        % transformacja cz�stotliwo�ci
        ws1t = (fs1*(fp2-fp1)) / (fs1^2 - fp1*fp2); 
        ws2t = (fs2*(fp2-fp1)) / (fs2^2 - fp1*fp2);
        ws = min(abs(ws1t), abs(ws2t));
    end
    
    % Przeliczenie na warto�� bezwzgl�dn�
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
        vc = 2*pi*sqrt(fp1*fp2);	% pulsacja �rodka
        % szeroko�� filtra wok� vc
        f_ps = [fp1,fp2,fs1,fs2];
        dv = 2*pi*(fp2-fp1);
        wzm_ps = [wzm_p, wzm_p, wzm_s, wzm_s];
        wzmdB_ps = [-apass, -apass, -astop, -astop];
    end
    
    disp([num2str(i) ' - FILTR: ' filtr ' - ' typ]); 
    
    % Obliczenie parametr�w N i w0 - funckja buttord
    wp = 1;
    N = ceil(log10( (10^(astop/10)-1) /(10^(apass/10)-1) ) / (2*log10(ws/wp)) );
    w0 = ws / (10^(astop/10)-1)^(1/(2*N));

    % Obliczenie biegun�w trans. prototypu - funcja buttap i zp2tf
    dfi0 = (2*pi)/(2*N);                    % k�t �kawa�ka tortu�
    fi = pi/2 + dfi0/2 + (0 : N-1)*dfi0;    % k�ty biegun�w
    p = w0*exp(j*fi);                       % bieguny
    z = [];                                 % zera
    wzm = real( prod(-p) );                 % wzmocnienie
    a = poly(p);         % bieguny --> wsp wielomianu mianownika A(z)
    b = wzm;             % wielomian licznika B(z)
%     z, p, b, a
    figure(4*i-3)
    plot( real(p), imag(p), 'x' ); grid; axis equal;
    title(['Po�o�enie biegun�w dla filtra:' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
    
    % Por�wnanie z Matlabem
    [NN,ww0] = buttord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystyk� cz�stotliwo�ciow� H(w)=B(w)/A(w)
    % zakres pulsacji unorm.; pulsacja granicy pasma przepustowego = 1
    w = 0 : 0.005 : 2;	
    H = freqs(b,a,w);   % alternatywa: H=polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i-2); subplot(211);
    plot(w,abs(H));
    grid; title(['Modu� prototypu LowPass dla' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]');
    
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(['Modu� prototypu LowPass w dB dla ' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    % Transformata cz�stotliwo�ci filtra analogowego: prototyp unormowany
    % --> wynikowy filtr
    
    % LowPass to LowPass: s=s/w0
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end
    % LowPass to HighPass: s=w0/s
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end
    % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end
    % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end
    
    b=wzm*poly(z); a=poly(p);
    % Poka� zera i bieguny po transformacji cz�stoliwo�ci
    figure(4*i-1)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(['Po�o�enie biegun�w dla filtra ' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
%     p, z
%     a, b
    printsys(b,a,'s');
    % Ko�cowa charakterystyka cz�stoliwo�ciowa
    NF = 1000;	% ile punkt�w
    fmin = 0;	% dolna cz�stotliwo��
    fmax = 5000;	% g�rna cz�stotliwo��
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie cz�stotliwo�ci
    w = 2*pi*f;     % wszystkie pulasacje
    H = freqs(b,a,w);	% alternatywa: H=polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i); subplot(211);
    plot( f,abs(H), f_ps, wzm_ps,'ro');
    grid; title(['Modu� dla filtra: ' filtr ' - ' typ]);
    xlabel('Czestotliwo�� [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro');
    axis([fmin,fmax,-100,20]);
    grid; title(['Modu� w dB filtra ' filtr ' - ' typ]);
    xlabel('Czestotliwo�� [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwo�� [Hz]'); ylabel('[rad]');
end

%% Filtry na podstawie prototypu Czebyszewa I typu
typ = "Czebyszewa I";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = 'LowPass';
        fpass = 1000; fstop = 4000; % Hz
        ws = fstop/fpass;   % czestotliwo�� znorm. s=s'/w0, w0=2*pi*fpass
    
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = 'HighPass';
        fstop = 2000; % cz�s. pasma przepustowego odpowiadaj�ca astop
        fpass = 3000; % cz�s. pasma zaporowego odpowiadaj�ca apass
        ws = fpass/fstop; % transformacja cz�s.i: s=w0/s', w0=2*pi*fpass
    
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = 'BandPass';
        fs1 = 1500; % dolna cz�stotliwo�� stop
        fp1 = 2000; % dolna cz�stotliwo�� pass
        fp2 = 3000; % g�rna cz�stotliwo�� pass
        fs2 = 3500; % g�rna cz�stotliwo�� stop
        % transformacja cz�stotliwo�ci
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    
    % filtr srodkowozaporowy BS
    else
        filtr = 'BandStop';
        fp1 = 1500; % dolna cz�stotliwo�� filtra pasmowego
        fs1 = 2000; % dolna cz�stotliwo�� filtra pasmowego
        fs2 = 3000; % g�rna cz�stotliwo�� filtra pasmowego
        fp2 = 3500; % g�rna cz�stotliwo�� filtra pasmowego
        % transformacja cz�stotliwo�ci
        ws1t = (fs1*(fp2-fp1)) / (fs1^2 - fp1*fp2); 
        ws2t = (fs2*(fp2-fp1)) / (fs2^2 - fp1*fp2);
        ws = min(abs(ws1t), abs(ws2t));
    end
    
    % Przeliczenie na warto�� bezwzgl�dn�
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
        vc = 2*pi*sqrt(fp1*fp2);	% pulsacja �rodka
        % szeroko�� filtra wok� vc
        f_ps = [fp1,fp2,fs1,fs2];
        dv = 2*pi*(fp2-fp1);
        wzm_ps = [wzm_p, wzm_p, wzm_s, wzm_s];
        wzmdB_ps = [-apass, -apass, -astop, -astop];
    end
    
    disp([num2str(i) ' - FILTR: ' filtr ' - ' typ]); 
    
    % Obliczenie parametr�w N i w0 
    wp = 1;
    Nreal = acosh(sqrt((10^(astop/10)-1) /...
        (10^(apass/10)-1))) / acosh(ws/wp);
    N = ceil(Nreal);
    epsi = sqrt(10^(apass/10)-1);
    D = asinh(1/epsi)/N; R1 = sinh(D); R2 = cosh(D);

    % Obliczenie biegun�w trans. prototypu - funcja zp2tf
    dfi0 = (2*pi)/(2*N);                    % k�t �kawa�ka tortu�
    fi = pi/2 + dfi0/2 + (0 : N-1)*dfi0;    % k�ty biegun�w
    p1 = R1 * exp(1i*fi);                    % bieguny R1
    p2 = R2 * exp(1i*fi);                    % bieguny R2
    p = real(p1) +1i*imag(p2);              % Polaczone bieguny
    z = [];                                 % zera
    wzm = prod(-p);                 % wzmocnienie
    a = poly(p);	% bieguny --> wsp wielomianu mianownika A(z)
    b = wzm;        % wielomian licznika B(z)
    if (rem(N,2)==0) b = b*10^(-apass/20); end

    figure(4*i-3+16)
    plot( real(p), imag(p), 'x' ); grid;
    title(['Po�o�enie biegun�w dla filtra:' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
    
    % Por�wnanie z Matlabem
    [NN,ww0] = cheb1ord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystyk� cz�stotliwo�ciow� H(w)=B(w)/A(w)
    % zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
    w = 0 : 0.005 : 2;
    H = freqs(b,a,w);	% alternatywa: H=polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i-2+16); subplot(211);
    plot(w,abs(H)); grid;
    title(['Modu� prototypu LowPass dla' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]');
    
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(['Modu� prototypu LowPass w dB dla ' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    % Transformata cz�stotliwo�ci filtra analogowego: prototyp unormowany
    % --> wynikowy filtr
    
    % LowPass to LowPass: s=s/w0
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end
    % LowPass to HighPass: s=w0/s
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end
    % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end
    % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end
    
    b=wzm*poly(z); a=poly(p);
    % Poka� zera i bieguny po transformacji cz�stoliwo�ci
    figure(4*i-1+16)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(['Po�o�enie biegun�w dla filtra ' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
%     p, z
%     a, b
    printsys(b,a,'s');
    % Ko�cowa charakterystyka cz�stoliwo�ciowa
    NF = 1000; % ile punkt�w
    fmin = 0; % dolna cz�stotliwo��
    fmax = 5000; % g�rna cz�stotliwo��
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie cz�stotliwo�ci
    w = 2*pi*f; % wszystkie pulasacje
    H = freqs(b,a,w); % alternatywa: H = polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i+16); subplot(211);
    plot(f, abs(H), f_ps, wzm_ps,'ro');
    grid; title(['Modu� dla filtra: ' filtr ' - ' typ]);
    xlabel('Czestotliwo�� [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro');
    axis([fmin,fmax,-100,20]);
    grid; title(['Modu� w dB filtra ' filtr ' - ' typ]);
    xlabel('Czestotliwo�� [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwo�� [Hz]'); ylabel('[rad]');
end

%% Filtry na podstawie prototypu Czebyszewa II typu
typ = "Czebyszewa II";
for i=1:4
    % filtr dolnoprzepustowy LP
    if i==1
        filtr = 'LowPass';
        fpass = 1000; fstop = 4000; % Hz
        ws = fstop/fpass;   % czestotliwo�� znorm. s=s'/w0, w0=2*pi*fpass
    
    % filtr gornoprzepustowy HP
    elseif i==2
        filtr = 'HighPass';
        fstop = 2000; % cz�s. pasma przepustowego odpowiadaj�ca astop
        fpass = 3000; % cz�s. pasma zaporowego odpowiadaj�ca apass
        ws = fpass/fstop; % transformacja cz�s.i: s=w0/s', w0=2*pi*fpass
    
    % filtr srodkowoprzepustowy BP
    elseif i==3
        filtr = 'BandPass';
        fs1 = 1500; % dolna cz�stotliwo�� stop
        fp1 = 2000; % dolna cz�stotliwo�� pass
        fp2 = 3000; % g�rna cz�stotliwo�� pass
        fs2 = 3500; % g�rna cz�stotliwo�� stop
        % transformacja cz�stotliwo�ci
        ws1t = (fs1^2 - fp1*fp2) / (fs1*(fp2-fp1));
        ws2t = (fs2^2 - fp1*fp2) / (fs2*(fp2-fp1));
        ws = min(abs(ws1t), abs(ws2t));
    
    % filtr srodkowozaporowy BS
    else
        filtr = 'BandStop';
        fp1 = 1500; % dolna cz�stotliwo�� filtra pasmowego
        fs1 = 2000; % dolna cz�stotliwo�� filtra pasmowego
        fs2 = 3000; % g�rna cz�stotliwo�� filtra pasmowego
        fp2 = 3500; % g�rna cz�stotliwo�� filtra pasmowego
        % transformacja cz�stotliwo�ci
        ws1t = (fs1*(fp2-fp1)) / (fs1^2 - fp1*fp2); 
        ws2t = (fs2*(fp2-fp1)) / (fs2^2 - fp1*fp2);
        ws = min(abs(ws1t), abs(ws2t));
    end
    
    % Przeliczenie na warto�� bezwzgl�dn�
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
        vc = 2*pi*sqrt(fp1*fp2);	% pulsacja �rodka
        % szeroko�� filtra wok� vc
        f_ps = [fp1,fp2,fs1,fs2];
        dv = 2*pi*(fp2-fp1);
        wzm_ps = [wzm_p, wzm_p, wzm_s, wzm_s];
        wzmdB_ps = [-apass, -apass, -astop, -astop];
    end
    
    disp([num2str(i) ' - FILTR: ' filtr ' - ' typ]); 
    
    % Obliczenie parametr�w N i w0
    wp = 1;
    Nreal = acosh(sqrt((10^(astop/10)-1) /...
        (10^(apass/10)-1))) / acosh(ws/wp);
    N = ceil(Nreal);
    epsi = sqrt(1 / (10^(apass/10)-1));
    D = asinh(1/epsi)/N; R1 = sinh(D); R2 = cosh(D);

    % Obliczenie biegun�w trans. prototypu - funcja buttap i zp2tf
    dfi0 = (2*pi)/(2*N);                    % k�t �kawa�ka tortu�
    fi = pi/2 + dfi0/2 + (0 : N-1)*dfi0;    % k�ty biegun�w
    p1 = R1 * exp(j*fi);                    % bieguny R1
    p2 = R2 * exp(j*fi);                    % bieguny R2
    p = real(p1) +1i*imag(p2);              % Polaczone bieguny
    z = 1i*sin(fi)                          % zera
    wzm = prod(-z) / prod(-p);              % wzmocnienie
    z = 1./z; p = 1./p;
    a = poly(p);        % bieguny --> wsp wielomianu mianownika A(z)
    b = wzm*poly(z);	% wielomian licznika B(z)

    figure(4*i-3+32)
    plot( real(p), imag(p), 'x' ); grid;
    title(['Po�o�enie biegun�w dla filtra:' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
    
    % Por�wnanie z Matlabem
    [NN,ww0] = cheb2ord( vp, vs, apass, astop, 's' );
	blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);
    
    % Oblicz charakterystyk� cz�stotliwo�ciow� H(w)=B(w)/A(w)
    % zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
    w = 0 : 0.005 : 2;
    H = freqs(b,a,w);	% alternatywa: H=polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i-2+32); subplot(211);
    plot(w,abs(H)); grid;
    title(['Modu� prototypu LowPass dla' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]');
    
    subplot(212); plot(w,20*log10(abs(H))); grid;
    title(['Modu� prototypu LowPass w dB dla ' filtr ' - ' typ]);
    xlabel('Pulsacja [rad/sek]'); ylabel('dB');
    % Transformata cz�stotliwo�ci filtra analogowego: prototyp unormowany
    % --> wynikowy filtr
    % LowPass to LowPass: s=s/w0
    if (i==1) [z,p,wzm] = lp2lpTZ(z,p,wzm,vp); end
    % LowPass to HighPass: s=w0/s
    if (i==2) [z,p,wzm] = lp2hpTZ(z,p,wzm,vp); end
    % LowPass to BandPass: s=(s^2+wc^2)/(dw*s)
    if (i==3) [z,p,wzm] = lp2bpTZ(z,p,wzm,vc,dv); end
    % LowPass to BandStop: s=(dw*s)/(s^2+wc^2)
    if (i==4) [z,p,wzm] = lp2bsTZ(z,p,wzm,vc,dv); end
    
    b=wzm*poly(z); a=poly(p);
    % Poka� zera i bieguny po transformacji cz�stoliwo�ci
    figure(4*i-1+32)
    plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
    title(['Po�o�enie biegun�w dla filtra ' filtr ' - ' typ]);
    xlabel('real'); ylabel('imag');
%     p, z
%     a, b
    printsys(b,a,'s');
    % Ko�cowa charakterystyka cz�stoliwo�ciowa
    NF = 1000; % ile punkt�w
    fmin = 0; % dolna cz�stotliwo��
    fmax = 5000; % g�rna cz�stotliwo��
    f = fmin : (fmax-fmin)/(NF-1) : fmax; % wszystkie cz�stotliwo�ci
    w = 2*pi*f; % wszystkie pulasacje
    H = freqs(b,a,w); % alternatywa: H=polyval( b,j*w)./polyval(a,j*w);
    
    figure(4*i+32); subplot(211);
    plot(f, abs(H), f_ps, wzm_ps,'ro');
    grid; title(['Modu� dla filtra: ' filtr ' - ' typ]);
    xlabel('Czestotliwo�� [Hz]');
    subplot(212);
    plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro');
    axis([fmin,fmax,-100,20]);
    grid; title(['Modu� w dB filtra ' filtr ' - ' typ]);
    xlabel('Czestotliwo�� [Hz]'); ylabel('dB');
    plot(f,unwrap(angle(H))); grid; title('Faza');
    xlabel('Czestotliwo�� [Hz]'); ylabel('[rad]');
end

%% Zaprojektowanie uk�adu elektronicznego dolnoprzepustowego filtra
% Butterwortha
% Dane projektowe
fpass = 8000;       % cz�sto. pasma przepustowego odpowiadaj�ca apass
fstop = 22050;      % cz�sto. pasma zaporowego odpowiadaj�ca astop
apass = 2;          % nieliniowo�� pasma przepustowego w dB (�zwis�)
astop = 40;         % t�umienie w pa�mie zaporowym

wzm_p = 10^(-apass/20);	% t�umienie pass -> wzmocnienie pass
wzm_s = 10^(-astop/20);	% t�umienie stop -> wzmocnienie stop
ws = fstop/fpass;   % transformacja cz�stotliwo�ci: s=s'/w0, w0=2*pi*fpass
vp = 2*pi*fpass; vs = 2*pi*fstop;
f_ps = [fpass, fstop];
wzm_ps = [wzm_p, wzm_s];
wzmdB_ps = [-apass, -astop];

wp = 1;
Nreal = log10( (10^(astop/10)-1)/(10^(apass/10)-1) ) / (2*log10(ws/wp));
N = ceil(Nreal);
w0 = ws / (10^(astop/10)-1)^(1/(2*N));

dfi0 = (2*pi)/(2*N);                    % k�t �kawa�ka tortu�
fi = pi/2 + dfi0/2 + (0 : N-1)*dfi0;	% k�ty biegun�w
p = w0*exp(1i*fi);                      % bieguny
z = [];                                 % zera
wzm = real(prod(-p));                   % wzmocnienie

figure(49)
plot( real(p), imag(p), 'x' ); grid;
title(['Po�o�enie biegun�w dla filtra dolnoprzepustowego']);
xlabel('real'); ylabel('imag');

b = wzm;                % wielomian licznika B(z)
a = poly(p);            % bieguny --> wsp wielomianu mianownika A(z)
printsys(b,a,'s');

% Por�wnanie z Matlabem
[NN,ww0] = buttord( vp, vs, apass, astop, 's' );
blad_N = N-NN; disp(['Blad rzedu: ' num2str(blad_N)]);

% Oblicz charakterystyk� cz�stotliwo�ciow� H(w)=B(w)/A(w)
% zakres pulsacji unormowanej; pulsacja granicy pasma przepustowego = 1
w = 0 : 0.005 : 2;
H = freqs(b,a,w);	% alternatywa: H=polyval( b,j*w)./polyval(a,j*w);

figure(50); subplot(211);
plot(w,abs(H)); grid; title(['Modu� prototypu LowPass']);
xlabel('Pulsacja [rad/sek]');

subplot(212); plot(w,20*log10(abs(H))); grid;
title(['Modu� prototypu LowPass w dB']);
xlabel('Pulsacja [rad/sek]'); ylabel('dB');

% Transformata cz�stotliwo�ci filtra analogowego: prototyp unormowany -->
% wynikowy filtr
[z,p,wzm] = lp2lpTZ(z,p,wzm,vp);     % LowPass to LowPass: s=s/w0

b=wzm*poly(z); a=poly(p);
% Poka� zera i bieguny po transformacji cz�stoliwo�ci
figure(51)
plot( real(z), imag(z), 'o',real(p),imag(p),'x' ); grid;
title(['Po�o�enie biegun�w dla filtra LowPass']);
xlabel('real'); ylabel('imag');
printsys(b,a,'s');

% Ko�cowa charakterystyka cz�stoliwo�ciowa
NF = 1000;	% ile punkt�w
fmin = 0;	% dolna cz�stotliwo��
fmax = 50000;	% g�rna cz�stotliwo��
f = fmin : (fmax-fmin)/(NF-1) : fmax;	% wszystkie cz�stotliwo�ci
w = 2*pi*f;     % wszystkie pulasacje
H = freqs(b,a,w);   % alternatywa: H=polyval( b,j*w)./polyval(a,j*w);

figure(52); subplot(211);
plot( f,abs(H), f_ps, wzm_ps,'ro');
grid; title(['Modu� dla filtra lowpass']);
xlabel('Czestotliwo�� [Hz]');
subplot(212);
plot(f,20*log10(abs(H)), f_ps, wzmdB_ps,'ro');
axis([fmin,fmax,-100,20]);
grid; title(['Modu� w dB filtra lowpass']);
xlabel('Czestotliwo�� [Hz]'); ylabel('dB');
plot(f,unwrap(angle(H))); grid; title('Faza');
xlabel('Czestotliwo�� [Hz]'); ylabel('[rad]');

% Oblicz elementy uk�adu ze wzmacniaczami operacyjnymi
p1 = [ p(1) conj(p(1)) ];
p2 = [ p(2) conj(p(2)) ];
p3 = p(4);
aw1 = poly(p1), aw2 = poly(p2), aw3 = poly(p3);
C = 10^(-9); RA=10^4; Rwy = 10^4;

disp('=== Uk�ad 1===')
a = aw1;
a2=a(1); a1=a(2); a0=a(3);
R = 1/(C*sqrt(a0))
RB = (2-a1/sqrt(a0)) * RA
K1 = 1+RB/RA

disp('=== Uk�ad 2 ===')
a = aw2;
a2=a(1); a1=a(2); a0=a(3);
R = 1/(C*sqrt(a0))
RB = (2-a1/sqrt(a0)) * RA
K2 = 1+RB/RA

disp('=== Uk�ad 3 ===')
a = aw3;
a1=a(1); a0=a(2);
R=1/(C*a0)
K3=1

disp('=== Obci��enie ===')
K=K1*K2*K3
G=1
Rx = (K/G)*Rwy
Ry = (G/K)/(1-G/K)*Rx



%%%%%%%%%%%%% DEFINICJE FUNKCJI ZAGNIEZDZONYCH %%%%%%%%%%%%%%

% Funkcja transformuj�ca filtr LP znromalizowany na wymagany LP
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

% Funkcja transformuj�ca filtr LP znromalizowany na wymagany HP
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


% Funkcja transformuj�ca filtr LP znromalizowany na wymagany BP
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


% Funkcja transformuj�ca filtr LP znromalizowany na wymagany BS
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