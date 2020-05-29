function x = sbp(w, A, t, n)
% Funcja generuj¹ca fale prostkatna bipolarna
% w - czêstoœæ, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t));
for i=1:length(n)
    for j=1:2:n(i)
        x(i,:) = x(i,:) + ((1/j)*sin(j*w*t));
    end
end
x = x*4*A/pi;
end


