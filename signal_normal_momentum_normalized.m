function value = signal_normal_momentum_normalized(signal, m)
%	signal_normal_momentum_normalized
%   Return normalized signal normal momentum
%   singal - signal, m - order
value = signal_normal_momentum(signal, m) / sum(signal);
end

