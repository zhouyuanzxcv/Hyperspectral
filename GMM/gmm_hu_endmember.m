function [E] = gmm_hu_endmember(I,A,D,w_jk,mu_jk,sigma_jk)
%GMM_HU_ENDMEMBER Estimated endmembers for each pixel based on estimated
%GMM parameters and abundances
% Input: 
%   I - rows by cols by B image cube
%   A - estimated N by M abundance matrix
%   D - noise covariance matrix
% Output:
%   E - M by B by N matrix of endmembers for each pixel

% if nargin < 7
%     options = [];
% end

[rows,cols,B] = size(I);
[N,M] = size(A);
Y = reshape_hsi(I);
% D = parse_param(options, 'D', 0.001^2*eye(size(I,3)));

disp('Start estimating endmembers for each pixel.');
t_start = tic;

% reduce the dimensions
Y_ori = Y;
reduced_dim = 30;
[Y, mapping] = pca(Y_ori, reduced_dim);
[mu_jk,sigma_jk] = project2ortho(mu_jk,sigma_jk,mapping.mean,mapping.M);
B = reduced_dim;
I = reshape(Y, [rows, cols, B]);
D = mapping.M'*D*mapping.M;

% prepare calculation
sigma_jk_inv = sigma_jk;
for j = 1:M
    for k = 1:size(sigma_jk{j}, 3)
        sigma_jk_inv{j}(:,:,k) = inv(sigma_jk{j}(:,:,k));
    end
end

E = E_init(I,A,mu_jk); % M by B by N matrix

if 0 % also estimate D
    % Sometimes the estimated D is singular, causing the estimation of E
    % involves a matrix of bad condition number
    errs = eval_obj_fun(Y,A,E,D,w_jk,mu_jk,sigma_jk);
    for iter = 1:20
        E = optimize_wrt_E(Y, A, E, D, w_jk, mu_jk, sigma_jk, sigma_jk_inv);
        D = optimize_wrt_D(Y, A, E);
        errs(end+1) = eval_obj_fun(Y,A,E,D,w_jk,mu_jk,sigma_jk);
        if test_convergence(errs, 0.001, 1)
            break;
        end
    end
    D = mapping.M*D*mapping.M' + 1e-18*eye(size(mapping.M,1));
else % use input D
    E = optimize_wrt_E(Y, A, E, D, w_jk, mu_jk, sigma_jk, sigma_jk_inv);
end 

% restore from projection
[~,~,E] = restore_from_projection([],[],E,mapping.mean,mapping.M);

t_elapsed = toc(t_start);
disp(['Elapsed time for estimating all the endmembers is ',num2str(t_elapsed)]);

end

function D = optimize_wrt_D(Y, A, E)
[N,B] = size(Y);
Y1 = zeros(N,B);

for i = 1:N
    Y1(i,:) = A(i,:) * E(:,:,i);
end
Y1 = Y - Y1;
D = (Y1' * Y1) / N;
% The original D is diagonal. But the projected D may not be diagonal
% D = diag(diag(D));

end

function E = optimize_wrt_E(Y, A, E, D, w_jk, mu_jk, sigma_jk, sigma_jk_inv)
[N,B] = size(Y);
M = size(A,2);

D_inv = inv(D);

max_iter = 50;
convergence_t = 0.001;

C = zeros(M*B, M*B);

errs = eval_obj_fun(Y,A,E,D,w_jk,mu_jk,sigma_jk);

for iter = 1:max_iter
    t_start = tic;
    % E-step
    gamma_njk = calc_gamma_njk(E, w_jk, mu_jk, sigma_jk);
    
    % M-step
    [C_nj,d_n] = calc_C_d(gamma_njk, mu_jk, sigma_jk_inv);
    for i = 1:N
        ds = d_n(:,i,:);
        for j = 1:M
            inds = ((j-1)*B+1:j*B);
            C(inds,inds) = C_nj(:,:,i,j);
        end
        E_i = D_inv * Y(i,:)' * A(i,:);
        E_i = (kron(A(i,:)'*A(i,:), D_inv) + C) \ (E_i(:) + ds(:));
        E(:,:,i) = reshape(E_i, [B M])';
    end
    
    % test convergence
    errs(end+1) = eval_obj_fun(Y,A,E,D,w_jk,mu_jk,sigma_jk);
    
    t_elapsed = toc(t_start);
    disp(['EM iteration ',num2str(iter),' in estimating E lasts ',...
        num2str(t_elapsed)]);
    
    if test_convergence(errs, convergence_t, 3)
        break;
    end
end

if 0
    figure('name','energy decrease for endmember estimation');
    plot(errs);
end


end

function [C_nj,d_n] = calc_C_d(gamma_njk, mu_jk, sigma_jk_inv)
M = length(gamma_njk);
N = size(gamma_njk{1}, 1);
B = size(mu_jk{1}, 2);
C_nj = zeros(B,B,N,M);
d_n = zeros(B,N,M);

for j = 1:M
    K_j = size(sigma_jk_inv{j}, 3);
    C = reshape(sigma_jk_inv{j}, [B*B K_j]) * gamma_njk{j}';
    C_nj(:,:,:,j) = reshape(C, [B B N]);
    Mu = zeros(B,K_j);
    for k = 1:K_j
        Mu(:,k) = sigma_jk_inv{j}(:,:,k) * mu_jk{j}(k,:)';
    end
    d_n(:,:,j) = Mu * gamma_njk{j}';
end
end

function val = eval_obj_fun(Y,A,E,D,w_jk,mu_jk,sigma_jk)
[N,M] = size(A);
B = size(Y,2);
%
val_prior = zeros(N,M);
for j = 1:M
    val_prior(:,j) = calc_log_gmm(squeeze(E(j,:,:))', w_jk{j}, mu_jk{j}, sigma_jk{j});
end

diff = Y;
for j = 1:M
    diff = diff - squeeze(E(j,:,:))' .* repmat(A(:,j), [1 B]);
end
[V,D1] = eig(inv(D));
F = diag(sqrt(diag(D1))) * V';
val_lsq = 0.5 * sum(sum((F * diff').^2, 1));

val = val_lsq - sum(sum(val_prior));

end

function E = E_init(I,A,mu_jk)
M = length(mu_jk);
K = zeros(1,M);
for j = 1:M
    K(j) = size(mu_jk{j},1);
end
K_all = K2K_all(K);

mu_all = calc_mu_all(mu_jk, K_all);
Y = reshape_hsi(I);

[N,B] = size(Y);
[K1,M] = size(K_all);

recon_error = zeros(N,K1);
mu_all1 = zeros(M,B,K1);
for k = 1:K1
    mu_all1(:,:,k) = mu_all(:,:,k)';
    recon_error(:,k) = sum((Y - A*mu_all1(:,:,k)).^2, 2);
end

[~,inds] = min(recon_error,[],2);
E = zeros(M,B,N);
for i = 1:N
    E(:,:,i) = mu_all1(:,:,inds(i));
end

end

