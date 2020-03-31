function rms = signal_rms(signal, dt)
%	signal_rms
%	Summary of this function goes here
%   Detailed explanation goes here
rms = sqrt( signal_power(signal, dt) );
end
