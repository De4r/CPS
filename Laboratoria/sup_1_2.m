function x = sup_1_2(w, A, t, n)
% Funcja generuj¹ca fale prostkatna unipolarna o wypelnieniu 1/2
% w - czêstoœæ, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t));
for i=1:length(n)
    for j=1:4:n(i)
        x(i,:) = x(i,:) + ((1/j)*cos(j*w*t));
    end
    for j=3:4:n(i)
        x(i,:) = x(i,:) - ((1/j)*cos(j*w*t));
    end
end
x = x*2*A/pi + A/2;
end


