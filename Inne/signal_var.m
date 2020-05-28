function var_value = signal_var(signal)
%	signal_var
%	Returns signal variance
var_value = 0;
signal_mean_value = signal_mean(signal);
for i=1:length(signal)
    var_value = var_value + (signal(i) - signal_mean_value)^2;
end
var_value = var_value / ( length(signal) - 1);
end
