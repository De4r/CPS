function x = sup_wyp(f, A, t, n, tau)
% Funcja generuj¹ca fale prostkatna unipolarna o dowolnym wypelnieniu
% w - czestotliwosc, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t)); T = 1/f;
for i=1:length(n)
    for j=1:n(i)
        x(i,:) = x(i,:) + sin(pi*j*tau/T)*cos(2*j*pi*f*t)/(pi*j*tau/T);
    end       
end
x = A*tau/T + 2*A*tau*x/T;
end

