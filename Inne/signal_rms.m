function rms = signal_rms(signal, dt)
%	signal_rms
%	Returns signal root mean square
rms = sqrt( signal_power(signal, dt) );
end
