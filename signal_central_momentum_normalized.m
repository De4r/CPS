function value = signal_central_momentum_normalized(signal, m)
%	signal_central_momentum_normalized
%   Normalized signal central momentum
%   signal - signal, m - order
nx1 = signal_normal_momentum_normalized(signal, 1); value = 0;
for i=1:length(signal)
    value = value + (i - nx1)^m * signal(i);
end
value = value / sum(signal);
end



