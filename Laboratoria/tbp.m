function x = tbp(w, A, t, n)
% Funcja generuj¹ca fale trojkatna bipolarna
% w - czestotliwosc, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t));
for i=1:length(n)
    for j=1:4:n(i)
        x(i,:) = x(i,:) + ((1/j^2)*sin(j*w*t));
    end
    for j=3:4:n(i)
        x(i,:) = x(i,:) - ((1/j^2)*sin(j*w*t));
    end      
end
x = x*8*A/(pi^2);
end
