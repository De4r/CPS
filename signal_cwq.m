function cwq = signal_cwq(signal)
%   center weight of signal quadrature
%   Summary of this function goes here
%   Detailed explanation goes here
licznik = 0; mianownik = 0;
for i=1:length(signal)
   licznik = licznik + i * signal(i)^2;
   mianownik = mianownik + signal(i)^2;
end
cwq = licznik/mianownik;
end

