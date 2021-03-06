function x = tupp(w, A, t, n)
% Funcja generująca fale trojkatna unipolarna pilokształtna
% w - częstość, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t));
for i=1:length(n)
    for j=1:n(i)
        x(i,:) = x(i,:) - ((1/j)*sin(j*w*t));
    end   
end
x = x*A/pi + A/2;
end

