function A = calc_A_from_mus(Y, mu_jk)
%CALC_A_FROM_MUS Summary of this function goes here
%   Detailed explanation goes here
M = length(mu_jk);
K = zeros(1,M);
for j = 1:M
    K(j) = size(mu_jk{j},1);
end
K_all = K2K_all(K);

mu_all = calc_mu_all(mu_jk, K_all);
[N,B] = size(Y);
[K1,M] = size(K_all);
recon_error = zeros(N,K1);
A_all = zeros(N,M,K1);

for k = 1:K1
    R = mu_all(:,:,k)';
    A = Y*R'*inv(R*R'+eye(M)*1e-6);
    A = project_to_simplex(A);
    recon_error(:,k) = sum((Y - A*R).^2, 2);
    A_all(:,:,k) = A;
end

[~,inds] = min(recon_error,[],2);
A = zeros(N,M);
for i = 1:N
    A(i,:) = A_all(i,:,inds(i));
end
