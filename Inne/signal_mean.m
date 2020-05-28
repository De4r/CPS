function mean = signal_mean(signal)
%   signal_mean
%   Return mean value of signal
mean = 0; % set variable to 0
for i=1:length(signal)
   mean = mean + signal(i); 
end
mean = mean/length(signal);
end

