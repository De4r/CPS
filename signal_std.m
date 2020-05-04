function std_value = signal_std(signal)
%	signal_std 
%	Retruns standard deviation
std_value = 0;
signal_mean_value = signal_mean(signal);
for i=1:length(signal)
    std_value = std_value + (signal(i) - signal_mean_value)^2;
end
std_value = std_value / ( length(signal) - 1);
std_value = sqrt(std_value);
end