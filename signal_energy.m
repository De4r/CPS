function energy = signal_energy(signal, dt)
%	signal_energy 
%	Summary of this function goes here
%   Detailed explanation goes here
energy = 0;
for i=1:length(signal)
    energy = energy + (signal(i))^2*dt;
end

end

