function sigma = snr2sigma(X,SNR)
%SNR2SIGMA Summary of this function goes here
%   Detailed explanation goes here
[N,B] = size(X);
p = mean(sum(X.^2,2));
sigma2 = p/(B*10^(SNR/10));
sigma = sqrt(sigma2);

end

