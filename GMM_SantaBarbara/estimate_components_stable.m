function [mapping,components] = estimate_components_stable(endmembers, options)
%ESTIMATE_COMPONENTS_STABLE Summary of this function goes here
%   Detailed explanation goes here
options = insert_param_when_absent(options,'reduced_dim',10);
mapping = calc_projection_from_library(endmembers, options.reduced_dim);  

disp('Run model selection multiple times to obtain stable GMM components');

% repeat multiple times
repeat_times = 15;
for i = 1:repeat_times
    [K,w_jk,mu_jk,sigma_jk] = estimate_components(endmembers, mapping, options);
    Ks{i} = K;
    w_jks{i} = w_jk;
    mu_jks{i} = mu_jk;
    sigma_jks{i} = sigma_jk;
end

% find the most frequent combination
K1s = zeros(1,repeat_times);
for i = 1:repeat_times
    K1s(i) = prod(Ks{i});
end
K1_most = mode(K1s);
inds = find(K1s == K1_most);

% find the most frequent w_jk
K_most = [];
w_jk_most = {};
mu_jk_most = {};
sigma_jk_most = {};
for j = 1:length(endmembers)
    x = zeros(1,length(inds));
    for i = 1:length(inds)
        x(i) = prod(w_jks{inds(i)}{j});
    end
    f = ksdensity(x,x);
    [~,ind_j] = max(f);
    ind_j = inds(ind_j);
    K_most(j) = Ks{ind_j}(j);
    w_jk_most{j} = w_jks{ind_j}{j};
    mu_jk_most{j} = mu_jks{ind_j}{j};
    sigma_jk_most{j} = sigma_jks{ind_j}{j};
end

disp(['The most frequent number of combinations is ',num2str(prod(K_most))]);
% return this combination
components = [];
components.K = K_most;
components.w_jk = w_jk_most;
components.mu_jk = mu_jk_most;
components.sigma_jk = sigma_jk_most;

end

