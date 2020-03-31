function min_value = signal_min(signal)
%	signal_min 
%	Summary of this function goes here
%   Detailed explanation goes here
min_value = 0;
for i=1:length(signal)
    if signal(i) < min_value
        min_value = signal(i);
    end
end
end

