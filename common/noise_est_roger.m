function sigmas = noise_est_roger(I)
%NOISE_EST_ROGER Estimate the standard deviation of noise for each band.
if ndims(I) == 3
    [rows,cols,B] = size(I);
    Y = reshape(I, [rows*cols, B]);
else
    Y = I;
end

[N,B] = size(Y);
Y1 = Y - repmat(mean(Y,1),N,1);
K = (1/N)*(Y1'*Y1);
[D,E] = decompose(inv(K));
sigmas = diag(D).^(-1);


function [D,E] = decompose(K)
d = sqrt(diag(K));
D = diag(d);
E = diag(1./d) * K * diag(1./d);
% mdiff(D*E*D,K)
