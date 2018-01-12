function mu_all = calc_mu_all(mu_jk, K_all)
% mu_all is B by M by K1
[~,B] = size(mu_jk{1});
[K1,M] = size(K_all);

mu_all = zeros(B,M,K1);
for i = 1:K1
    for j = 1:M
        mu_all(:,j,i) = mu_jk{j}(K_all(i,j),:)';
    end
end


end

