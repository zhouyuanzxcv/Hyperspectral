function [ output_args ] = gmm_hu(I,M,options)
%GMM_HU Summary of this function goes here
%   Detailed explanation goes here
beta1 = 1e4;
beta2 = 0; % 2e4
eta = 0.05;

rho1 = 1e4;

if isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end

names = cell(1,M);
for i = 1:M
    names{i} = ['endmember ',num2str(i)];
end

%% initialize parameters
[Y,~,rows,cols] = reshape_hsi(I,[]);
[N,B] = size(Y);

beta1 = beta1*B/M;
beta2 = beta2*B/M;
rho1 = rho1*N/M^2;

if beta1 < 1e-9
    beta1 = 1e-9;
end

[W,Neighbors] = image2graph(I,eta,1e-9);

D = diag(sum(W,2));
L = D - W;
L = sparse(L);
KL = L - beta2/beta1*speye(N);

W = ones(M,M);
D = diag(sum(W,2));
H = D - W;

[mu_jk,sigma_jk,w_jk,K,A] = gmm_init(I,M);

%% Create all the k indices
% K = [1 2 3 1];
K_inds  = cell(1,M);
for j = 1:M
    K_inds{j} = (1:K(j));
end
K_all = cartprod(K_inds{:});
K1 = size(K_all,1);


%% iterate by MM (EM)
w_k = w_jk2w_k(w_jk,K_all);

sigma = Y_noise;

I_B = eye(B);

delta_t0 = 1e-12;

delta_t_mu = delta_t0;
delta_t_sigma = delta_t0;
delta_t_A = delta_t0;

der_mu0 = mu_jk;
der_sigma0 = sigma_jk;
der_A0 = zeros(N,M);
for j = 1:length(der_mu0)
    der_mu0{j}(:,:) = 0;
    der_sigma0{j}(:,:,:) = 0;
end

s = [];
s.N = N;
s.B = B;
s.K1 = K1;
s.sigma = sigma;
s.K_all = K_all;
s.Y = Y;
s.beta1 = beta1;
s.KL = KL;

s.A = A;
s.mu_jk = mu_jk;
s.sigma_jk = sigma_jk;
s.w_k = w_k;

Is = (1:N*B);
Is = repmat(reshape(Is, [B,1,N]), 1, B);
Is = Is(:);

Js = (1:N*B);
Js = repmat(reshape(Js, [1,B,N]), B, 1);
Js = Js(:);

s.Is = Is;
s.Js = Js;

max_iter = 200;
eval_totals = zeros(max_iter, 1);
eval_Ns = zeros(max_iter, 1);
eval_As = zeros(max_iter, 1);
eval_Rs = zeros(max_iter, 1);

for iter = 1:200    
    %% E step
    % update gamma_nk
    
    N_nk = calc_gaussians(A, sigma, mu_jk, sigma_jk, K_all, Y, s);
    
    gamma_nk = (ones(N,1)*w_k) .* N_nk;
    gamma_nk = gamma_nk ./ repmat(sum(gamma_nk,2), 1, K1);
    s.gamma_nk = gamma_nk;
    
    %% M step
    % update w_k
    w_k = sum(gamma_nk, 1) / N;
    s.w_k = w_k;
    
    delta_t = delta_t0;
    disp('Process M step optimization.');
    

    % update mu
    der_mu = calc_der_mu(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y);
    s.der_mu = der_mu;
    delta_t_mu = calc_time_step_adaptive(@eval_obj_fun_mu, @update_mu, ...
        mu_jk, s, delta_t_mu, delta_t0);
    mu_jk = update_mu(mu_jk, s, delta_t_mu);
    s.mu_jk = mu_jk;

    % update sigma
    der_sigma = calc_der_sigma(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y);
    s.der_sigma = der_sigma;
    delta_t_sigma = calc_time_step_adaptive(@eval_obj_fun_sigma, @update_sigma, ...
        sigma_jk, s, delta_t_sigma, delta_t0);
    sigma_jk = update_sigma(sigma_jk, s, delta_t_sigma);
    s.sigma_jk = sigma_jk;

    % update A
    der_A = calc_der_A(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y, KL, beta1);
    s.der_A = der_A;
    delta_t_A = calc_time_step_adaptive(@eval_obj_fun_A, @update_A, ...
        A, s, delta_t_A, delta_t0);
    A = update_A(A, s, delta_t_A);
    s.A = A;
    
    if mod(iter,40) == 0
        endmember_scatter_plot_end_var(Y,w_jk,mu_jk,sigma_jk,names);
        set(gcf,'name',['Scatter plot of iteration ',num2str(iter)]);
        show_abundances(A,size(I,1),size(I,2));
        set(gcf,'name',['Abundances of iteration ',num2str(iter)]);
        pause(0.01);
    end
    
    [eval_total, eval_N, eval_A, eval_R] = calc_obj_fun(A, ...
        sigma, mu_jk, sigma_jk, K_all, Y, KL, beta1, w_k, s);
    eval_totals(iter) = eval_total;
    eval_Ns(iter) = eval_N;
    eval_As(iter) = eval_A;
    eval_Rs(iter) = eval_R;
    
    disp(['EM iteration ', num2str(iter)]);

end

w_jk = w_k2w_jk(w_k, K_all);

show_abundances(A,size(I,1),size(I,2));
figure('name', 'Total objective function value vs iteration number');
plot(eval_totals);

figure('name', 'Data fidelity term value vs iteration number');
plot(eval_Ns);

figure('name', 'Abundance smoothness term value vs iteration number');
plot(eval_As);

figure('name', 'Endmembers clossness term value vs iteration number');
plot(eval_Rs);


function [eval_total, eval_N, eval_A, eval_R] = calc_obj_fun(A, ...
    sigma, mu_jk, sigma_jk, K_all, Y, KL, beta1, w_k, options)
N = size(A,1);

N_nk = calc_gaussians(A, sigma, mu_jk, sigma_jk, K_all, Y, options);
eval_N = -sum(log(sum((ones(N,1)*w_k) .* N_nk, 2)));
eval_A = (beta1/2)*trace(A'*KL*A);
eval_R = 0;
eval_total = eval_N + eval_A + eval_R;


function N_nk = calc_gaussians(A, sigma, mu_jk, sigma_jk, K_all, Y, options)
[N,B] = size(Y);
K1 = size(K_all,1);

[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);
% old implementation

% N_nk1 = zeros(N,K1);
% for n = 1:N
%     for k = 1:K1
%         y_n_mu_nk = Y(n,:)'-mu_nk(:,:,n,k);
%         N_nk1(n,k) = det(sigma_nk(:,:,n,k))^(-1/2) * ...
%             exp(-0.5 * y_n_mu_nk' * (sigma_nk(:,:,n,k) \ y_n_mu_nk));
%     end
% end
% N_nk1 = N_nk1 / (2*pi)^(B/2);
    
% new implementation

N_nk = zeros(N,K1);
for k = 1:K1
    y = logmvn(Y, mu_nk(:,:,k), sigma_nk(:,:,:,k), options);
    N_nk(:,k) = exp(y);
end 

% norm(N_nk - N_nk1) / norm(N_nk)


function val = eval_obj_fun_mu(params, options)
mu_jk = params;

A = options.A;
sigma = options.sigma;
% mu_jk = options.mu_jk;
sigma_jk = options.sigma_jk;
K_all = options.K_all;
gamma_nk = options.gamma_nk;
Y = options.Y;
rho1 = options.rho1;
w_k = options.w_k;
% beta1 = options.beta1;
% KL = options.KL;
H = options.H;

[N,B] = size(Y);
K1 = size(K_all,1);
[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);


val2 = zeros(1,K1);
for k = 1:K1
    Y1 = (Y - mu_nk(:,:,k))';
    Y2 = Y1 .* repmat(gamma_nk(:,k)', B, 1);
    sigma_k = block_diag(sigma_nk(:,:,:,k), options);
    val2(k) = Y1(:)' * (sigma_k \ Y2(:));
end
val = sum(sum(0.5 * val2));

% tic
% N_nk = zeros(N,K1);
% for k = 1:K1
%     N_nk(:,k) = logmvn(Y, reshape(mu_nk(:,:,:,k),B,N)', sigma_nk(:,:,:,k), options);
% end 
% 
% val3 = -sum(sum(0.5 * gamma_nk .* N_nk));
% toc

val_mu = 0;
[mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all);
for k = 1:K1
    val_mu = val_mu + w_k(k) * trace(mu_all(:,:,k) * H * mu_all(:,:,k)');
end
val_mu = rho1 / 2 * val_mu;

val = val + val_mu;


function val = eval_obj_fun_sigma(params, s)
sigma_jk = params;

A = s.A;
sigma = s.sigma;
mu_jk = s.mu_jk;
% sigma_jk = options.sigma_jk;
K_all = s.K_all;
gamma_nk = s.gamma_nk;
Y = s.Y;
% beta1 = options.beta1;
% KL = options.KL;

N = size(A,1);
K1 = size(K_all,1);

[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);
% val1 = zeros(N,K1);
% for n = 1:N
%     for k = 1:K1
%         y_n_mu_nk = Y(n,:)' - mu_nk(:,:,n,k);
%         val1(n,k) = log(det(sigma_nk(:,:,n,k))) + y_n_mu_nk' * (sigma_nk(:,:,n,k) \ y_n_mu_nk);
%     end
% end
% 
% val = sum(sum(0.5 * gamma_nk .* val1));

N_nk = zeros(N,K1);
for k = 1:K1
    N_nk(:,k) = logmvn(Y, mu_nk(:,:,k), sigma_nk(:,:,:,k), s);
end 

val = -sum(sum(0.5 * gamma_nk .* N_nk));



function val = eval_obj_fun_A(params, options)
A = params;

% A = options.A;
sigma = options.sigma;
mu_jk = options.mu_jk;
sigma_jk = options.sigma_jk;
K_all = options.K_all;
gamma_nk = options.gamma_nk;
Y = options.Y;
beta1 = options.beta1;
KL = options.KL;

N = size(A,1);
K1 = size(K_all,1);
[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);
% val1 = zeros(N,K1);
% for n = 1:N
%     for k = 1:K1
%         y_n_mu_nk = Y(n,:)' - mu_nk(:,:,n,k);
%         val1(n,k) = log(det(sigma_nk(:,:,n,k))) + y_n_mu_nk' * (sigma_nk(:,:,n,k) \ y_n_mu_nk);
%     end
% end

% val = sum(sum(0.5 * gamma_nk .* val1)) + beta1/2 * trace(A'*KL*A);

N_nk = zeros(N,K1);
for k = 1:K1
    N_nk(:,k) = logmvn(Y, mu_nk(:,:,k), sigma_nk(:,:,:,k), options);
end 

val = -sum(sum(0.5 * gamma_nk .* N_nk)) + beta1/2 * trace(A'*KL*A);




function der_mu = calc_der_mu(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y)
[N,M] = size(A);
B = size(mu_jk{1},2);
K1 = size(K_all,1);
K = max(K_all,[],1);

[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);
lambda_nk = calc_lambda_nk(gamma_nk,mu_nk,sigma_nk,Y);

% calculate der_mu
% tic
% der_mu = mu_jk;
% for j = 1:M
%     for k = 1:K(j)
%         lambda_sum = sum(lambda_nk(:,:,K_all(:,j)==k), 3);
%         lambda_alpha = repmat(A(:,j), [1,B]);
%         der_mu{j}(k,:) = -sum(lambda_sum .* lambda_alpha, 1);
%     end
% end
% toc

% tic
der_mu = mu_jk;
temp = zeros(M,B,K1);
for k = 1:K1
    temp(:,:,k) = A' * lambda_nk(:,:,k);
end

for j = 1:M
    for k = 1:K(j)
        der_mu{j}(k,:) = -sum(temp(j,:,K_all(:,j)==k), 3)';
    end
end
% toc

% mdif(der_mu{1},der_mu1{1});

function der_sigma = calc_der_sigma(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y)
[N,M] = size(A);
[~,B] = size(Y);
K = max(K_all,[],1);
K1 = size(K_all,1);

[mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, sigma, mu_jk, sigma_jk, K_all);
[lambda_nk,psi_nk] = calc_lambda_psi_nk(gamma_nk,mu_nk,sigma_nk,prec,Y);

% calculate der_mu, der_sigma, der_A
% tic
% der_sigma1 = sigma_jk;
% for j = 1:M
%     for k = 1:K(j)
%         psi_sum = sum(psi_nk(:,:,:,K_all(:,j)==k), 4);
%         psi_alpha = reshape(repmat(A(:,j)'.^2, [B*B,1]), [B,B,N]);
%         der_sigma1{j}(:,:,k) = -sum(psi_sum .* psi_alpha, 3);        
%     end
% end
% toc

% tic
psi_k = reshape(psi_nk, [B*B,N,K1]);
der_sigma = sigma_jk;
temp = zeros(M,B*B,K1);
for k = 1:K1
    temp(:,:,k) = (A.^2)' * psi_k(:,:,k)';
end

for j = 1:M
    for k = 1:K(j)
        der_sigma{j}(:,:,k) = reshape(-sum(temp(j,:,K_all(:,j)==k), 3)', B, B);
    end
end
% toc

% disp('');

function lambda_nk = calc_lambda_nk(gamma_nk,mu_nk,sigma_nk,Y)
[N,B,K1] = size(mu_nk);

% tic
% lambda_nk1 = zeros(N,B,K1);
% for n = 1:N
%     for k = 1:K1
%         lambda_nk1(n,:,k) = gamma_nk(n,k) * ...
%             (sigma_nk(:,:,n,k) \ (Y(n,:)' - mu_nk(n,:,k)'));
%     end
% end
% toc

% tic
lambda_nk = zeros(N,B,K1);
for k = 1:K1
    Y1 = Y - mu_nk(:,:,k);
    Y1 = Y1';
    Y2 = Y1 .* repmat(gamma_nk(:,k)', B, 1);
    sigma_k = block_diag(sigma_nk(:,:,:,k));
    lambda_nk(:,:,k) = reshape(sigma_k \ Y2(:), [B,N])';
end
% toc

function [lambda_nk,psi_nk] = calc_lambda_psi_nk(A, sigma, mu_jk, ...
    sigma_jk, K_all, gamma_nk, Y)
[mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, sigma, mu_jk, sigma_jk, K_all);

[N,B,K1] = size(mu_nk);

% tic
% lambda_nk1 = zeros(N,B,K1);
% psi_nk1 = zeros(B,B,N,K1);
% 
% for k = 1:K1
%     Y1 = Y - mu_nk(:,:,k);
%     for n = 1:N
%         sigma_y_n_mu_nk = prec(:,:,n,k) * Y1(n,:)';
%         lambda_nk1(n,:,k) = gamma_nk(n,k) * sigma_y_n_mu_nk;
%         psi_nk1(:,:,n,k) = 0.5 * gamma_nk(n,k) * ( ...
%             sigma_y_n_mu_nk * sigma_y_n_mu_nk' - prec(:,:,n,k));
%     end
% end
% toc

% tic
lambda_nk = zeros(N,B,K1);
psi_nk = zeros(B,B,N,K1);

for k = 1:K1
    Y1 = Y - mu_nk(:,:,k);
    sigma_y_n_mu_nk = multiprod(prec(:,:,:,k), Y1', [1 2], [1]);
    lambda_nk(:,:,k) = (repmat(gamma_nk(:,k)', B, 1) .* sigma_y_n_mu_nk)';
    
    tmp1 = reshape(sigma_y_n_mu_nk, [B 1 N]);
    tmp2 = reshape(sigma_y_n_mu_nk, [1 B N]);
    tmp3 = multiprod(tmp1, tmp2, [1 2], [1 2]) - prec(:,:,:,k);
    psi_nk(:,:,:,k) = 0.5 * multiprod(reshape(gamma_nk(:,k), [1 1 N]), tmp3);
end
% toc
% mdiff(lambda_nk,lambda_nk1);
% mdiff(psi_nk,psi_nk1);


function der_A = calc_der_A(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, ...
    Y, KL, beta1)
[N,K1] = size(gamma_nk);
M = length(mu_jk);
B = size(mu_jk{1},2);

[mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, sigma, mu_jk, sigma_jk, K_all);
[lambda_nk,psi_nk] = calc_lambda_psi_nk(gamma_nk,mu_nk,sigma_nk,prec,Y);


% calculate der_mu, der_sigma, der_A
% tic
% 
% der_A1 = zeros(N,M);
% for n = 1:N
%     for j = 1:M
%         for k = 1:K1
%             der_A1(n,j) = der_A1(n,j) - sum(lambda_nk(n,:,k) .* mu_jk{j}(K_all(k,j),:)) ...
%                 - 2 * A(n,j) * sum(sum(psi_nk(:,:,n,k) .* sigma_jk{j}(:,:,K_all(k,j))));
%         end
%     end
% end
% toc

% tic
[mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all);
psi_k = reshape(psi_nk, [B*B,N,K1]);

der_A = zeros(N,M);
temp1 = zeros(N,M);
for k = 1:K1
    der_A = der_A - lambda_nk(:,:,k) * mu_all(:,:,k);
    temp1 = temp1 - psi_k(:,:,k)' * sigma_all(:,:,k);
end
der_A = der_A + 2 * A .* temp1;
% toc

der_A = der_A + beta1 * KL * A;


function mu_jk_new = update_mu(mu_jk, options, delta_t)
der_mu = options.der_mu;

mu_jk_new = mu_jk;

for j = 1:length(mu_jk)
    mu_jk_new{j} = mu_jk{j} - delta_t * der_mu{j};
end

function sigma_jk_new = update_sigma(sigma_jk, options, delta_t)
der_sigma = options.der_sigma;

sigma_jk_new = sigma_jk;

for j = 1:length(sigma_jk)
    sigma_jk_new{j} = sigma_jk{j} - delta_t * der_sigma{j};
%     sigma_jk_new{j} = sigma_jk{j};
    for k = 1:size(sigma_jk{j},3)
        cov_mat = sigma_jk_new{j}(:,:,k);
        [V,D] = eig((cov_mat + cov_mat')/2);
        d = diag(D);
        if (min(d) < 0)
            disp('There is a negative eigenvalue in the updated sigma_jk');
            d(d<0) = 1e-6 * max(abs(d));
            D = diag(d);
        end
        sigma_jk_new{j}(:,:,k) = V*D*V';
    end
end

function A_new = update_A(A, options, delta_t)
der_A = options.der_A;

A_new = A - delta_t * der_A;
A_new = project_to_simplex(A_new);

function [mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all)
[N,M] = size(A);
K1 = size(K_all,1);
B = size(mu_jk{1},2);

[mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all);

sigma2_I = sigma^2 * eye(B);

% old implementation
% mu_nk1 = zeros(B,1,N,K1);
% sigma_nk1 = zeros(B,B,N,K1);
% for n = 1:N
%     for k = 1:K1
%         mu_nk1(:,:,n,k) = mu_all(:,:,k) * A(n,:)';
%         sigma_nk1(:,:,n,k) = reshape(sigma_all(:,:,k) * (A(n,:).^2)', B, B) + sigma2_I;
%     end
% end

% new implementation is 10 times faster
mu_nk = zeros(N,B,K1);
sigma_nk = zeros(B*B,N,K1);
for k = 1:K1
    mu_nk(:,:,k) = A * mu_all(:,:,k)';
    sigma_nk(:,:,k) = sigma_all(:,:,k) * (A.^2)' + repmat(sigma2_I(:),1,N);
end
% mu_nk = reshape(mu_nk,[B,1,N,K1]); 
sigma_nk = reshape(sigma_nk,[B,B,N,K1]);

function [mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all)
[~,B] = size(mu_jk{1});
[K1,M] = size(K_all);

mu_all = zeros(B,M,K1);
sigma_all = zeros(B*B,M,K1);
for i = 1:K1
    for j = 1:M
        mu_all(:,j,i) = mu_jk{j}(K_all(i,j),:)';
        sigma_all(:,j,i) = reshape(sigma_jk{j}(:,:,K_all(i,j)), B*B, 1);
    end
end


function [mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, sigma, mu_jk, sigma_jk, K_all)
[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);
% [N,B,K1] = size(mu_nk);
% tic
% prec1 = zeros(size(sigma_nk));
% for n = 1:N
%     for k = 1:K1
%         prec1(:,:,n,k) = inv(sigma_nk(:,:,n,k));
%     end
% end
% toc
% 
% tic
prec = multinv(sigma_nk);
% toc
% 
% mdiff(prec,prec1);

% tic
% prec1 = cell(1,K1);
% for k = 1:K1
%     sigma_k = block_diag(sigma_nk(:,:,:,k), options);
%     R = chol(sigma_k);
%     S = inv(R);
%     prec1{k} = S*S';
% end
% 
% prec2 = zeros(B,B,N,K1);
% for k = 1:K1
%     [~,~,S] = find(prec1{k});
%     prec2(:,:,:,k) = reshape(S, [B,B,N]);
% end
% 
% toc
% 
% mdif(prec,prec2)


function w_k = w_jk2w_k(w_jk, K_all)
w_k = ones(1, size(K_all,1));
for i = 1:size(K_all,1)
    k = K_all(i,:);
    for j = 1:length(w_jk)
        w_k(i) = w_k(i)*w_jk{j}(k(j));
    end
end

function w_jk = w_k2w_jk(w_k, K_all)
M = size(K_all,2);
w_jk = cell(1,M);
K = zeros(1,M);
for j = 1:M
    K(j) = max(K_all(:,j));
    w_jk{j} = zeros(1,K(j));
end

for j = 1:M
    for l = 1:K(j)
        w = w_k(K_all(:,j)==l);
        w_jk{j}(l) = sum(w);
    end
end

%% obsolete
function delta_t = calc_time_step(A, sigma, mu_jk, sigma_jk, K_all, ...
    gamma_nk, Y, delta_t, delta_t0, der_mu, der_sigma, der_A, beta1, KL)
val_ori = eval_obj_fun_M(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y, beta1, KL);
[mu_jk_new, sigma_jk_new, A_new] = update_params(mu_jk, sigma_jk, A, der_mu, der_sigma, der_A, delta_t);
val_new = eval_obj_fun_M(A_new, sigma, mu_jk_new, sigma_jk_new, K_all, gamma_nk, Y, beta1, KL);

if val_new < val_ori
    val_old = val_new;
    while 1
        delta_t = delta_t * 10;
        [mu_jk_new, sigma_jk_new, A_new] = update_params(mu_jk, sigma_jk, A, ...
            der_mu, der_sigma, der_A, delta_t);
        val = eval_obj_fun_M(A_new, sigma, mu_jk_new, sigma_jk_new, K_all, gamma_nk, Y, beta1, KL);
        if val < val_old
            val_old = val;
        else
            delta_t = delta_t / 10;
            break;
        end
    end
else
    while delta_t > delta_t0
        delta_t = delta_t / 10;
        [mu_jk_new, sigma_jk_new, A_new] = update_params(mu_jk, sigma_jk, A, ...
            der_mu, der_sigma, der_A, delta_t);
        val = eval_obj_fun_M(A_new, sigma, mu_jk_new, sigma_jk_new, K_all, gamma_nk, Y, beta1, KL);
        if val < val_ori
            break;
        end
    end
end


function [mu_jk_new, sigma_jk_new, A_new] = update_params(mu_jk, sigma_jk, A, ...
    der_mu, der_sigma, der_A, delta_t)
mu_jk_new = mu_jk;
sigma_jk_new = sigma_jk;

for j = 1:length(mu_jk)
    mu_jk_new{j} = mu_jk{j} - delta_t * der_mu{j};
    sigma_jk_new{j} = sigma_jk{j} - delta_t * der_sigma{j};
%     sigma_jk_new{j} = sigma_jk{j};
    for k = 1:size(sigma_jk{j},3)
        cov_mat = sigma_jk_new{j}(:,:,k);
        [V,D] = eig((cov_mat + cov_mat')/2);
        d = diag(D);
        if (min(d) < 0)
            disp('There is a negative eigenvalue in the updated sigma_jk');
            d(d<0) = 1e-6 * max(abs(d));
            D = diag(d);
        end
        sigma_jk_new{j}(:,:,k) = V*D*V';
    end
end

A_new = A - delta_t * der_A;
A_new = project_to_simplex(A_new);


function val = eval_obj_fun_M(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y, beta1, KL)
N = size(A,1);
K1 = size(K_all,1);
[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);
val1 = zeros(N,K1);
for n = 1:N
    for k = 1:K1
        y_n_mu_nk = Y(n,:)' - mu_nk(:,:,n,k);
        val1(n,k) = log(det(sigma_nk(:,:,n,k))) + y_n_mu_nk' * (sigma_nk(:,:,n,k) \ y_n_mu_nk);
    end
end

val = sum(sum(0.5 * gamma_nk .* val1)) + beta1/2 * trace(A'*KL*A);


function [der_mu, der_sigma, der_A] = calc_derivatives(A, sigma, mu_jk, ...
    sigma_jk, K_all, gamma_nk, Y, beta1, KL)
[N,M] = size(A);
B = size(mu_jk{1},2);
K1 = size(K_all,1);
K = max(K_all,[],1);

der_mu = mu_jk;
der_sigma = sigma_jk;

[mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, sigma, mu_jk, sigma_jk, K_all);

lambda_nk = zeros(B,1,N,K1);
psi_nk = zeros(B,B,N,K1);
der_A = zeros(N,M);
for n = 1:N
    for k = 1:K1
        sigma_y_n_mu_nk = prec(:,:,n,k) * (Y(n,:)' - mu_nk(:,:,n,k));
        lambda_nk(:,:,n,k) = gamma_nk(n,k) * sigma_y_n_mu_nk;
        psi_nk(:,:,n,k) = 0.5 * gamma_nk(n,k) * ( ...
            sigma_y_n_mu_nk * sigma_y_n_mu_nk' - prec(:,:,n,k));
    end
end
% calculate der_mu, der_sigma, der_A
for j = 1:M
    for k = 1:K(j)
        lambda_sum = sum(lambda_nk(:,:,:,K_all(:,j)==k), 4);
        lambda_alpha = reshape(repmat(A(:,j)', [B,1]), [B,1,N]);
        der_mu{j}(k,:) = -sum(lambda_sum .* lambda_alpha, 3);

        psi_sum = sum(psi_nk(:,:,:,K_all(:,j)==k), 4);
        psi_alpha = reshape(repmat(A(:,j)'.^2, [B*B,1]), [B,B,N]);
        der_sigma{j}(:,:,k) = -sum(psi_sum .* psi_alpha, 3);        
    end
end

for n = 1:N
    for j = 1:M
        for k = 1:K1
            der_A(n,j) = der_A(n,j) - sum(lambda_nk(:,:,n,k) .* mu_jk{j}(K_all(k,j),:)') ...
                - 2 * A(n,j) * sum(sum(psi_nk(:,:,n,k) .* sigma_jk{j}(:,:,K_all(k,j))));
        end
    end
end
der_A = der_A + beta1 * KL * A;


function [der_mu, der_sigma] = calc_derivative_sigma(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y)
[N,M] = size(A);
B = size(mu_jk{1},2);
K1 = size(K_all,1);
K = max(K_all,[],1);

der_mu = mu_jk;
der_sigma = sigma_jk;

for j = 1:M
    der_mu{j} = zeros(size(der_mu{j}));
end

[mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, sigma, mu_jk, sigma_jk, K_all);

psi_nk = zeros(B,B,N,K1);
for n = 1:N
    for k = 1:K1
        sigma_y_n_mu_nk = prec(:,:,n,k) * (Y(n,:)' - mu_nk(:,:,n,k));
        psi_nk(:,:,n,k) = 0.5 * gamma_nk(n,k) * ( ...
            sigma_y_n_mu_nk * sigma_y_n_mu_nk' - prec(:,:,n,k));
    end
end
% calculate der_sigma
for j = 1:M
    for k = 1:K(j)
        psi_sum = sum(psi_nk(:,:,:,K_all(:,j)==k), 4);
        psi_alpha = reshape(repmat(A(:,j)'.^2, [B*B,1]), [B,B,N]);
        der_sigma{j}(:,:,k) = -sum(psi_sum .* psi_alpha, 3);
        
%         temp = zeros(B,B);
%         for n = 1:N
%             for k1 = 1:K1
%                 if k == K_all(k1,j)
%                     temp = temp - A(n,j)^2 * psi_nk(:,:,n,k1);
%                 end
%             end
%         end
    end
end
