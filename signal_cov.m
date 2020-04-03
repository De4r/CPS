function c = signal_cov(x, y, scale)
%   signal_cov
%   Summary of this function goes here
%   Detailed explanation goes here
%   x - sygnal x, y - sygnal y
%   scale 0: - brak
%         1: - normalizacja obciazona
%         2: - normalizacja nieobciazona
c = zeros(length(x)-1);
c = signal_corr(x-mean(x), y-mean(y), scale);
end

