function [f, M, W] = fft_from_signal(y, fs)
%	fft_from_signal 
%   Summary of this function goes here
%   Detailed explanation goes here
%   y - signal matrix, with signals as rows
%   f - frequency, M - power sepctrum, W - spectrum
N = length(y);
for i=1:size(y, 1)
    fft_moc=fft(y(i, 1:N));
    moc_wid=fft_moc.*conj(fft_moc)/N;
    widmo=sqrt(fft_moc.*conj(fft_moc))/N;
    f=fs*(0:N/2-1)/N;
    M(i,:)=moc_wid;
    W(i,:)=widmo;
end
M = M(:, floor(1:N/2));
W = W(:, floor(1:N/2));
end

