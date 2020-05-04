function value = signal_central_momentum(signal, m)
%	signal_central_momentum
%   signal - signal, m - order

value = 0; kx1 = signal_normal_momentum(signal, 1);
for i=1:length(signal)
    value = value + (i - kx1)^m * signal(i);
end


