function power = signal_power(signal, dt)
%	signal_power 
%	Summary of this function goes here
%   Detailed explanation goes here
power = signal_energy(signal, dt)/length(signal);
end