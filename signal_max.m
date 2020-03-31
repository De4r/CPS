function max_value = signal_max(signal)
%	signal_max 
%	Summary of this function goes here
%   Detailed explanation goes here
max_value = 0;
for i=1:length(signal)
    if signal(i) > max_value
        max_value = signal(i);
    end
end
end

