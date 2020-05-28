function c = signal_cov(x, y, scale)
%   signal_cov
%   
%   Auto and cross covariance if passed signals
%   x - signal x, y - signal y (sygnaly)
%   scale 0: - none (brak)
%         1: - biased (normalizacja obciazona)
%         2: - unbiased (normalizacja nieobciazona)
c = zeros(length(x)-1);
c = signal_corr(x-mean(x), y-mean(y), scale);
end

