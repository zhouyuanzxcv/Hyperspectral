function K_all = K2K_all(K)
%% Transform between K and K_all

M = length(K);
K_inds  = cell(1,M);
for j = 1:M
    K_inds{j} = (1:K(j));
end
K_all = cartprod(K_inds{:});

