function x = swd(w, A, t, n)
% Funcja generuj¹ca fale sinusoidalna wyprostowana dwupo³owkow¹
% w - czestotliwosc, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t));
for i=1:length(n)
    for j=1:n(i)
        x(i,:) = x(i,:) + (1/(4*j^2-1))*cos(2*j*w*t);
    end   
end
x = x*(-4)*A/pi + 2*A/pi;
end

