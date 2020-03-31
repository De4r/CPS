function vsqacwq = signal_vsqacwq(signal)
%   variance of signal quadrature around center of quadarture signal weight
%   xD
%   Summary of this function goes here
%   Detailed explanation goes here
cwq = signal_cwq(signal);
licznik = 0; mianownik = 0;
for i=1:length(signal)
   licznik = licznik + (i - cwq)^2 * signal(i)^2;
   mianownik = mianownik + signal(i)^2;
end
vsqacwq = licznik / mianownik;
end

