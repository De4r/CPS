function x = swj(w, A, t, n)
% Funcja generuj¹ca fale sinusoidalna wyprostowana jednopo³owkow¹
% w - czêstoœæ, A - amplituda
% t - wektor czasu, n - rzad ciagu
x=zeros(length(n), length(t));
for i=1:length(n)
    for j=1:n(i)
        x(i,:) = x(i,:) + (1/(4*j^2-1))*cos(2*j*w*t);
    end
    x(i,:) = x(i,:)*A*(-2)/pi + A/pi + sin(w*t)*A/2;
end
end

