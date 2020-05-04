function value = signal_normal_momentum(signal, m)
%	signal_normal_momentum
%   Return signal normal momentum
%   singal - signal, m - order
value = 0;
for i=1:length(signal)
    value = value + i^m*signal(i);
end

