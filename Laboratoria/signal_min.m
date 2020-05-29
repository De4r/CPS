function min_value = signal_min(signal)
%	signal_min 
%	Return min value of signal
min_value = 0;
for i=1:length(signal)
    if signal(i) < min_value
        min_value = signal(i);
    end
end
end

