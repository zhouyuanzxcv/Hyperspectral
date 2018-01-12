function K = K_all2K(K_all)
%% Transform between K and K_all

M = size(K_all,2);
K = zeros(1,M);
for j = 1:M
    K(j) = max(K_all(:,j));
end