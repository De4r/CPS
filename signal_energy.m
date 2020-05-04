function energy = signal_energy(signal, dt)
%	Signal energy
%	signal - signal, dt - time step
energy = 0;
if nargin > 1
  dt = dt;
else
  dt = 1;
end
for i=1:length(signal)
    energy = energy + (signal(i))^2*dt;
end

end

