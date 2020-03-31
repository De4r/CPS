function value = signal_normal_momentum_normalized(signal, m)
%	signal_normal_momentum_normalized
%   Summary of this function goes here
%   Detailed explanation goes here
value = signal_normal_momentum(signal, m) / sum(signal);
end

