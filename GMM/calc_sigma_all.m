function sigma_all = calc_sigma_all(sigma_jk, K_all)
% sigma_all is B^2 by M by K1
B = size(sigma_jk{1},1);
[K1,M] = size(K_all);

sigma_all = zeros(B*B,M,K1);
for i = 1:K1
    for j = 1:M
        sigma_all(:,j,i) = reshape(sigma_jk{j}(:,:,K_all(i,j)), B*B, 1);
    end
end


end

