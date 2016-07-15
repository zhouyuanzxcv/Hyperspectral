function y = logmvn(X, Mu, Sigma, options)
%LOGMVN Evaluate the logarithms of normal density functions.
%The function evaluates the logarithm of Gaussian density function given 
%data at X(n,:) with center at Mu(n,:) and covariance matrix Sigma(:,:,n).
%The code can handle different means and covariance matrices for each 
%sample. It is much faster than the native function mvnpdf.
%
% Input: 
%   X - N by D data. N is the number of samples. D is the number of
%       dimensions.
%   Mu - N by D means.
%   Sigma - D by D by N covariance matrices.
%   Options - structure containing indices for transforming the Sigma to a
%       sparse matrix. Usually no need to set.
% Output:
%   y - N by 1 results that contain the logarithm of the Gaussians.
%
% Example:
%
% N = 1000;
% D = 2;
% Mu = [1 -1]; 
% Mu = repmat(Mu, 1000, 1) + randn(N,D);
% SN = repmat(reshape(abs(randn(N,1)), [1,1,N]), D, D) .* repmat(eye(D), [1,1,N]);
% Sigma = [.9 .4; .4 .3];
% Sigma = repmat(Sigma, [1,1,N]) + SN;
% X = mvnrnd(Mu,Sigma,1000); 
% 
% tic
% y = logmvn(X, Mu, Sigma);
% y = exp(y);
% toc
% 
% tic
% y1 = mvnpdf(X, Mu, Sigma);
% toc
% 
% norm(y-y1)
%
% Author: Yuan Zhou (zhouyuanzxcv@gmail.com)

% Copyright (C) 2015


[N,B] = size(X);

%% Transform the sigma matrix to a sparse block diagonal matrix
if nargin > 3
    Is = options.Is;
    Js = options.Js;
else
    Is = (1:N*B);
    Is = repmat(reshape(Is, [B,1,N]), 1, B);
    Is = Is(:);

    Js = (1:N*B);
    Js = repmat(reshape(Js, [1,B,N]), B, 1);
    Js = Js(:);
end

% Construct the sparse matrix
Sigma1 = sparse(Is,Js,reshape(Sigma,N*B*B,1),N*B,N*B);

%% Compute the Cholesky decomposition
Y1 = X - Mu;
[R,err] = chol(Sigma1);
if err ~= 0
    error('A Covariance matrix is not positive definite.');
end

x1 = reshape(Y1', 1, N*B) / R;
x1 = sum(reshape(x1, B, N).^2, 1)';
x2 = sum(reshape(log(diag(R)),B,N), 1)';

y = -0.5 * x1 - x2 - B * log(2*pi) / 2;

end

