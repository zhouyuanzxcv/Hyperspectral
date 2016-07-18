function sigmas = noise_est_roger(Y)
%NOISE_EST_ROGER Estimate the standard deviation of noise for each band.
if ndims(Y) == 3
    [rows,cols,B] = size(Y);
    Y = reshape(Y, [rows*cols, B]);
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
