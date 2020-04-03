function x = signal_ifft(X)
%	signal_ifft
%   Summary of this function goes here
%   Detailed explanation goes here
%   X -signal, x - signal in time
%   For length N input vector x, the DFT is a length N vector X,
%   with elements
%                    N
%      X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
%                   n=1
%   The inverse DFT (computed by IFFT) is given by
%                    N
%      x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N.
%                   k=1

x = zeros(1, length(X));
for n=1:length(x)
    temp = 0;
    for k=1:length(x)
        temp = temp + X(k)*exp(1i*2*pi*(k-1)*(n-1)/length(x));
    end
    x(1, n) = temp/length(x);
end