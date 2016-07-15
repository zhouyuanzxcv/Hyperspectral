function SNR = calc_snr(Y,A,R)
%CALC_SNR Summary of this function goes here
%   Detailed explanation goes here
A = double(A);
X = A*R;
W = Y-X;
SNR = mean(sum(X.^2,2)) / mean(sum(W.^2,2));
SNR = 10*log10(SNR);

end

