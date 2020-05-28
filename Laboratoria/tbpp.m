function x = tbpp(w, A, t, n)
% Funcja generuj¹ca fale trojkatna bipolarna pilokszta³tna
% w - czestotliwosc, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t));
for i=1:length(n)
    for j=1:2:n(i)
        x(i,:) = x(i,:) + ((1/j)*sin(j*w*t));
    end
    for j=2:2:n(i)
        x(i,:) = x(i,:) - ((1/j)*sin(j*w*t));
    end      
end
x = x*2*A/pi;
end

