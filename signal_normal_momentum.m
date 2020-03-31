function value = signal_normal_momentum(signal, m)
%	signal_normal_momentum
%   Summary of this function goes here
%   Detailed explanation goes here
value = 0;
for i=1:length(signal)
    value = value + i^m*signal(i);
end

