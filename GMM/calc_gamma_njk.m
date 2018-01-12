function gamma_njk = calc_gamma_njk(E, w_jk, mu_jk, sigma_jk)
[M,B,N] = size(E);
gamma_njk = cell(1,M);
for j = 1:M
    M_j = squeeze(E(j,:,:))';
    K_j = length(w_jk{j});
    N_jk = zeros(N,K_j);
    for k = 1:K_j
        N_jk(:,k) = logmvn(M_j, mu_jk{j}(k,:), sigma_jk{j}(:,:,k));
    end
    gamma_njk{j} = calc_gamma_E_step(N_jk, w_jk{j});
end

end