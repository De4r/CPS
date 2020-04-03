function r = signal_corr(x, y, scale)
%   signal_corr
%   Summary of this function goes here
%   Detailed explanation goes here
%   x - sygnal x, y - sygnal y
%   scale 0: - brak
%         1: - normalizacja obciazona
%         2: - normalizacja nieobciazona
if size(x) == size(y);
	N = length(x);
	r = zeros(size(x));
    r2 = zeros(size(x));
	if scale == 0;
        for k=0:N-1
            r(k+1) = sum( x(1:N-k).*conj(y(1+k:N)));
            r2(k+1) = sum( x(1+k:N).*conj(y(1:N-k)));
        end
    elseif scale == 1;
        for k=0:N-1
            r(k+1) = sum( x(1:N-k).*conj(y(1+k:N))) / N;
            r2(k+1) = sum( x(1+k:N).*conj(y(1:N-k))) / N;
        end
    else scale == 2;
        for k=0:N-1
            r(k+1) = sum( x(1:N-k).*conj(y(1+k:N))) / (N-k);
            r2(k+1) = sum( x(1+k:N).*conj(y(1:N-k))) / (N-k);
        end
	end
end
%r = [fliplr(r) r(2:end)];
r = [fliplr(r) r2(2:end)];
end