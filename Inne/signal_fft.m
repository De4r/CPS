function X = signal_fft(x)
%	signal_fft
%   FFT of signal
%   x -signal, xf - fourier result
%   For length N input vector x, the DFT is a length N vector X,
%   with elements
%                    N
%      X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
%                   n=1
%   The inverse DFT (computed by IFFT) is given by
%                    N
%      x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N.
%                   k=1
assert(isvector(x));

N = length(x);

for k=1:N
	temp = 0;
	for n=1:N
        temp = temp + x(n)*exp(-2i*pi*(k-1)*(n-1)/N);
    end
	X(k) = temp;
end

end

