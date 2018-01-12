function value = calc_log_gmm(X, w, Mu, Sigma, options)
%CALC_LOG_GMM Calculate the logarithm of the GMM.
% Input:
%   X: N by B data
%   w: 1 by K vector, weights of the Gaussians
%   Mu: K by B (single mean) or N by B by K (multiple means)
%   Sigma: B by B by K (single covariance matrix) by B by B by N by K
%
% Output:
%   value: logarithm of the GMM
if nargin < 5
    options = [];
end

K = length(w);
[N,B] = size(X);

if size(Mu,1) == size(X,1) && size(Mu,3) == K
    is_single = 0;
elseif size(Mu,1) == K
    is_single = 1;
end

N_nk = zeros(N,K);

if is_single
    for k = 1:K
%         N_nk(:,k) = log(mvnpdf(X, Mu(k,:), Sigma(:,:,k)));
        N_nk(:,k) = loggausspdf(X, Mu(k,:), Sigma(:,:,k));
    end
else
    for k = 1:K
        N_nk(:,k) = logmvn(X, Mu(:,:,k), Sigma(:,:,:,k), options);
    end
end

% Need to avoid too large negative logarithm
N_nk_max = max(N_nk,[],2);

% eval_N1 = -sum(log(sum((ones(N,1)*w_k) .* exp(N_nk), 2)));
eval_N = log(sum(repmat(w,N,1) .* exp(N_nk - repmat(N_nk_max,1,K)), 2));
value = N_nk_max + eval_N;


function y = loggausspdf(X, mu, Sigma)
X = X';
mu = mu';
d = size(X,1);
X = bsxfun(@minus,X,mu);
[U,p]= cholcov(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
Q = U'\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(U)));   % normalization constant
y = -(c+q)/2;

y = y';