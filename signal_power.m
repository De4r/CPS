function power = signal_power(signal, dt)
%	signal_power 
%	Return singal power 
power = signal_energy(signal, dt)/length(signal);
end