function [K,w_jk,mu_jk,sigma_jk] = estimate_components(endmembers, mapping, options)
%ESTIMATE_COMPONENTS Summary of this function goes here
%   Detailed explanation goes here
options = insert_param_when_absent(options, 'max_num_comp', 4);

M = length(endmembers);

A = [];
X = [];
for j = 1:M
    N1 = size(endmembers{j},1);
    X1 = gmm_project(endmembers{j}, mapping);
    A1 = zeros(N1,M);
    A1(:,j) = 1;
    A = cat(1,A,A1);
    X = cat(1,X,X1);
end

sizes = [0 0];
[K,w_jk,mu_jk,sigma_jk,A1] = estimate_num_comp(X, A, sizes, 0, options.max_num_comp, options);


end

