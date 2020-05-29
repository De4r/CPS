function x = tup(w, A, t, n)
% Funcja generuj¹ca fale trojkatna unipolarna
% w - czêstoœæ, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t));
for i=1:length(n)
    for j=0:n(i)
        x(i,:) = x(i,:) + ((1/((2*j+1)^2))*cos((2*j+1)*w*t));
    end    
end
x = x*(-4*A)/(pi^2) + A/2;
end

