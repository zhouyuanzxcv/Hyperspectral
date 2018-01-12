function gamma_nk = calc_gamma_E_step(N_nk, w_k)
%CALC_GAMMA Calculate gamma in E step while scaling both the denominator
% and numerator by the largest quantity in N.
% Input:
%   N_nk - logrithm of Gaussian evaluations (N by K)
%   w_k - prior probabilities (1 by K)
% Output:
%   gamma_nk - soft membership in E-step (N by K)
K1 = length(w_k);
N = size(N_nk, 1);
max_N_nk = max(N_nk, [], 2);
N_nk1 = exp(N_nk - repmat(max_N_nk, [1,K1])); % avoid too large negative 

gamma_nk = (ones(N,1)*w_k) .* N_nk1;
gamma_nk = gamma_nk ./ repmat(sum(gamma_nk,2), 1, K1);


end

