function [A,R,w_jk,mu_jk,sigma_jk,extra] = gmm_hu(I,M,options)
%GMM_HU Summary of this function goes here
%   Detailed explanation goes here
beta1 = 5;
beta2 = 5;
rho1 = 0;
rho2 = 0;

eta = 0.05;

A_gt = [];
fix_A = 0;
fix_w_k = 0;
fix_sigma_jk = 0;
fix_mu_jk = 0;

convergence_thresh = 0.002;

shrink_size = 2;
max_num_comp = 4;
beta2_decay = 0.05;

use_last_init = 0;

skip_phase3 = 0;
skip_phase2 = 1;
reduced_dim = 10;
max_iter = 200;
show_fig = 1;

use_mex_optimization = 1;


project_mapping = [];

D = 0.001^2 * eye(size(I,3));

names = cell(1,M);
for i = 1:M
    names{i} = ['endmember ',num2str(i)];
end

if isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end

if max(I(:)) > 100 % I could be transformed by PCA before input, hence > 1
    disp('Warning! I is not in the range 0 - 1');
end

start_t = tic;

%% Project to low dimensional space
I_ori = I;
[Y_ori,~,rows,cols] = reshape_hsi(I_ori,[]);
if isempty(project_mapping)
    [Y, mapping] = pca(Y_ori, reduced_dim);
else
    mapping = project_mapping;
%     reduced_dim = length(mapping.mean);
    Y = (Y_ori - repmat(mapping.mean, size(Y_ori,1), 1)) * mapping.M;
end
% Y = gmm_project(Y_ori, mapping);
[N,B] = size(Y);
I = reshape(Y, [rows, cols, B]);
sigma = mapping.M'*D*mapping.M;

%% initialize parameters

beta1 = beta1*B/M;
beta2 = beta2*B/M;
rho1 = rho1*N/M^2;
rho2 = rho2*N/(M*B);

if beta1 < 1e-9
    beta1 = 1e-9;
end

[L,KL,H] = calc_Laplacians(I, eta, beta1, beta2, M);

if use_last_init
    load('init.mat');
    for j = 1:M
        sigma_jk{j}(:,:,1) = sigma0^2 * eye(B);
    end
elseif fix_mu_jk && fix_sigma_jk % supervised unmixing
%     R = gmm2mean(w_jk, mu_jk);
%     A = Y*R'*inv(R*R'+eye(M)*1e-6);
%     A = project_to_simplex(A);
    A = calc_A_from_mus(Y, mu_jk);
    
    K = options.K;
    w_jk = options.w_jk;
    mu_jk = options.mu_jk;
    sigma_jk = options.sigma_jk;
else % unsupervised unmixing
    [mu_jk,sigma_jk,w_jk,K,A] = gmm_init(I,M,options);
    save('init.mat','mu_jk','sigma_jk','w_jk','K','A');
end

%% Create all the k indices
K_all = K2K_all(K);
K1 = size(K_all,1);

%% iterate by MM (EM)
w_k = w_jk2w_k(w_jk,K_all);

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

options = [];
options.sigma = sigma;
options.K_all = K_all;
options.Y = Y;
options.beta1 = beta1;
options.rho1 = rho1;
options.rho2 = rho2;
options.KL = KL;
options.H = H;

options.A = A;
options.mu_jk = mu_jk;
options.sigma_jk = sigma_jk;
options.w_k = w_k;

Is = (1:N*B);
Is = repmat(reshape(Is, [B,1,N]), 1, B);
Is = Is(:);

Js = (1:N*B);
Js = repmat(reshape(Js, [1,B,N]), B, 1);
Js = Js(:);

options.Is = Is;
options.Js = Js;
if use_mex_optimization
    options.sigma_update_in_logmvn = sparse(Is,Js,ones(N*B*B,1),N*B,N*B);
end

eval_totals = [];
eval_Ns = [];
eval_As = [];
eval_Rs = [];
eval_Ss = [];


if fix_A
    A_gt = double(A_gt);
    P = permute_abundances(A, A_gt);
    A = A_gt * P;
%     A = A_gt;
    options.A = A;
end

start_sigma = 0;
sigma_iter = Inf;

if fix_mu_jk && fix_sigma_jk
    start_split = 1;
    split_iter = 0;
else
    start_split = 0;
    split_iter = Inf;
end

beta1_ori = beta1;
beta2_ori = beta2;
rho1_ori = rho1;
rho2_ori = rho2;

test_conv_w = 5;

epc = EvolveProcessCapture();
epc = epc.begin(max_iter);

[beta1,beta2,rho1,rho2,options] = rescale_params(sigma_jk, beta1_ori, ...
    beta2_ori, rho1_ori, rho2_ori, options);

for iter = 1:max_iter    
%     [beta1,beta2,rho1,rho2,options] = rescale_params(sigma_jk, beta1_ori, ...
%     beta2_ori, rho1_ori, rho2_ori, options);

    %% E step
    % update gamma_nk
    tstart = tic;
    
    N_nk = calc_log_gaussians(A, sigma, mu_jk, sigma_jk, K_all, Y, options);
    gamma_nk = calc_gamma_E_step(N_nk, w_k);
    options.gamma_nk = gamma_nk;
    
    %% M step
    % update w_k
    if ~fix_w_k
        w_k = sum(gamma_nk, 1) / N;
        options.w_k = w_k;
    end
        
    % update mu
    if ~fix_mu_jk
        der_mu = calc_der_mu(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y, rho1, H, w_k);
        options.der_mu = der_mu;
        delta_t_mu = calc_time_step_adaptive(@eval_obj_fun_mu, @update_mu, ...
            mu_jk, options, delta_t_mu, delta_t0);
        mu_jk = update_mu(mu_jk, options, delta_t_mu);
        options.mu_jk = mu_jk;
    end

    % update sigma
    if ~fix_sigma_jk
        if ~start_sigma && test_convergence(eval_totals, convergence_thresh, ...
                test_conv_w)
            if ~skip_phase2
                start_sigma = 1;
                sigma_iter = iter;
                disp('Start estimating the covariance matrices.');
            else
                start_sigma = 1;
                sigma_iter = -Inf;
            end
        end
        if start_sigma && ~skip_phase2
            der_sigma = calc_der_sigma(A, sigma, mu_jk, sigma_jk, ...
                K_all, gamma_nk, Y, rho2, w_k);
            options.der_sigma = der_sigma;
            delta_t_sigma = calc_time_step_adaptive(@eval_obj_fun_sigma, @update_sigma, ...
                sigma_jk, options, delta_t_sigma, delta_t0);
        %     delta_t_sigma = delta_t_sigma * max(0, (1-1/(0.1*iter)));
            sigma_jk = update_sigma(sigma_jk, options, delta_t_sigma);
            options.sigma_jk = sigma_jk;
        end
    end

    % update A
    if ~fix_A
        der_A = calc_der_A(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y, KL, beta1);
        options.der_A = der_A;
        delta_t_A = calc_time_step_adaptive(@eval_obj_fun_A, @update_A, ...
            A, options, delta_t_A, delta_t0);
        A = update_A(A, options, delta_t_A);
        options.A = A;
    end
    
    %% Show intermediate results and test convergence    
    [eval_total, eval_N, eval_A, eval_R, eval_S] = calc_obj_fun(A, sigma, ...
        mu_jk, sigma_jk, K_all, Y, KL, beta1, H, rho1, rho2, w_k, options);
    eval_totals = [eval_totals; eval_total];
    eval_Ns = [eval_Ns; eval_N];
    eval_As = [eval_As; eval_A];
    eval_Rs = [eval_Rs; eval_R];
    eval_Ss = [eval_Ss; eval_S];
        
    if ~start_split && start_sigma && iter > sigma_iter + 2*test_conv_w && ...
            test_convergence(eval_totals, convergence_thresh, test_conv_w)
        start_split = 1;
%         start_sigma = 1;
        split_iter = iter;
        
        if ~skip_phase3
            disp('Start splitting the centers');

            [K,w_jk,mu_jk,sigma_jk,A1] = estimate_num_comp(Y, A, ...
                [rows,cols], shrink_size, max_num_comp);
            if show_fig
                show_intermediate_abundances(A1, rows, cols, iter);
                set(gcf,'name','Abundances for estimating the number of components');
            end

            K_all = K2K_all(K);
            w_k = w_jk2w_k(w_jk, K_all);
    %         w_k1 = sort(w_k,2,'descend');
    %         figure('name','w_k values'); plot(w_k1);

            K1 = size(K_all,1);
            options.K_all = K_all;
            options.w_k = w_k;
            options.mu_jk = mu_jk;
            options.sigma_jk = sigma_jk;

            if 0
                opts = struct('show_approx',1,'w_jk',{w_jk},'mu_jk',...
                    {mu_jk},'sigma_jk',{sigma_jk});
                if show_fig
                    hist_end_var(Y,A1,names,1,opts);
                end
            end

            KL = L - beta2_decay * beta2/beta1 * speye(N);
            options.KL = KL;

            fix_w_k = 1;
            fix_mu_jk = 1;
            fix_sigma_jk = 1;
        end
    end
    
    if show_fig
        epc = epc.updateEndmVar(Y,A,rows,cols,w_jk,mu_jk,sigma_jk,names,iter);
    end
    
    telapsed = toc(tstart);
    if mod(iter, 10) == 0
        disp(['EM iteration ', num2str(iter), ' lasts ', num2str(telapsed)]);
    end
    
    if iter > split_iter + test_conv_w*2 && test_convergence( ...
            eval_totals, convergence_thresh, test_conv_w)
        break;
    end
end

% if show_fig
%     show_intermediate_abundances(A, rows, cols, iter);
% end

w_jk = w_k2w_jk(w_k, K_all);

if show_fig
    figure('name', 'Energy value trends');
    subplot(2,3,1); plot(eval_totals); xlabel('Total obj fun vs iter number');
    subplot(2,3,2); plot(eval_Ns); xlabel('Data fidelity vs iter number');
    subplot(2,3,3); plot(eval_As); xlabel('Abundance smoothness vs iter number');
    subplot(2,3,4); plot(eval_Rs); xlabel('Endmembers clossness vs iter number');
    subplot(2,3,5); plot(eval_Ss); xlabel('Covariance matrix size vs iter number');
end

% Restore to original space
[mu_jk,sigma_jk] = gmm_restore_from_projection(mu_jk,sigma_jk,mapping);
R = gmm2mean(w_jk, mu_jk);

epc = epc.end(iter);

extra = [];
extra.eval_totals = eval_totals;
extra.eval_Ns = eval_Ns;
extra.eval_As = eval_As;
extra.eval_Rs = eval_Rs;
extra.eval_Ss = eval_Ss;
extra.frames_scatter = epc.FramesScatter;
extra.frames_abund = epc.FramesAbund;

elapsed_t = toc(start_t);
disp(['Total algorithm execution time is ',num2str(elapsed_t)]);


function [L,KL,H] = calc_Laplacians(I, eta, beta1, beta2, M)
[rows,cols,B] = size(I);
N = rows * cols;

[W,Neighbors] = image2graph(I,eta,1e-9);

D = diag(sum(W,2));
L = D - W;
L = sparse(L);
KL = L - beta2/beta1 * speye(N);

W = ones(M,M);
D = diag(sum(W,2));
H = D - W;


function [beta1,beta2,rho1,rho2,options] = rescale_params(sigma_jk, ...
    beta1, beta2, rho1, rho2, options)
max_eig = -Inf;
for j = 1:length(sigma_jk)
    for k = 1:size(sigma_jk{j},3)
        d = eig(sigma_jk{j}(:,:,k));
        max_eig = max(max(d),max_eig);
    end
end
beta1 = beta1 / max_eig;
beta2 = beta2 / max_eig;
rho1 = rho1 / max_eig;
rho2 = rho2 / max_eig;
options.beta1 = beta1;
options.beta2 = beta2;
options.rho1 = rho1;
options.rho2 = rho2;

function show_intermediate_abundances(A, rows, cols, iter)
show_abundances(A,rows,cols,['Abundances of iteration ',num2str(iter)]);
pause(0.01);
jf = get(handle(gcf),'JavaFrame');
jf.setMinimized(1);

% function conv = test_convergence(eval_totals, size, threshold)
% if length(eval_totals) <= size * 2
%     conv = 0;
%     return;
% end
% initial_slope = eval_totals(1) - eval_totals(2);
% prev = mean(eval_totals(end - size*2 + 1: end - size));
% curr = mean(eval_totals(end - size + 1: end));
% conv = (prev - curr) < threshold * initial_slope;

%% Evaluate the total objective function
function [eval_total, eval_N, eval_A, eval_R, eval_S] = calc_obj_fun(A, ...
    sigma, mu_jk, sigma_jk, K_all, Y, KL, beta1, H, rho1, rho2, w_k, options)
[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);
eval_N = -sum(calc_log_gmm(Y, w_k, mu_nk, sigma_nk, options));

eval_R = calc_prior_mu(w_k, mu_jk, K_all, H, rho1);
eval_S = calc_prior_sigma(w_k, sigma_jk, K_all, rho2);

eval_A = (beta1/2) * trace(A'*KL*A);

eval_total = eval_N + eval_A + eval_R + eval_S;

function val_mu = calc_prior_mu(w_k, mu_jk, K_all, H, rho1)
K1 = size(K_all,1);
val_mu = 0;
mu_all = calc_mu_all(mu_jk, K_all);
for k = 1:K1
    val_mu = val_mu + w_k(k) * trace(mu_all(:,:,k) * H * mu_all(:,:,k)');
end
val_mu = rho1 / 2 * val_mu;

function val_sigma = calc_prior_sigma(w_k, sigma_jk, K_all, rho2)
[~,M] = size(K_all);
val2 = 0;
K = K_all2K(K_all);
w_jk = w_k2w_jk(w_k,K_all);
for j = 1:M
    for k = 1:K(j)
        val2 = val2 + w_jk{j}(k) * trace(sigma_jk{j}(:,:,k)'*sigma_jk{j}(:,:,k));
%         val2 = val2 + trace(sigma_jk{j}(:,:,k)'*sigma_jk{j}(:,:,k));
    end
end

% sigma_all = calc_sigma_all(sigma_jk, K_all);
% K1 = size(K_all,1);
% val3 = zeros(1,K1);
% for k = 1:K1
%     val3(k) = trace(sigma_all(:,:,k)' * sigma_all(:,:,k));
% end
% val3 = sum(w_k .* val3);
% mdiff(val2,val3)

val_sigma = rho2/2 * val2;

function N_nk = calc_log_gaussians(A, sigma, mu_jk, sigma_jk, K_all, Y, options)
[N,B] = size(Y);
K1 = size(K_all,1);

[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);

N_nk = zeros(N,K1);
for k = 1:K1
    N_nk(:,k) = logmvn(Y, mu_nk(:,:,k), sigma_nk(:,:,:,k), options);
end 

%% Evaluate the objective functions for projected gradient descent
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

val_mu = calc_prior_mu(w_k, mu_jk, K_all, H, rho1);

val = val + val_mu;


function val = eval_obj_fun_sigma(params, options)
sigma_jk = params;

A = options.A;
sigma = options.sigma;
mu_jk = options.mu_jk;
% sigma_jk = options.sigma_jk;
K_all = options.K_all;
gamma_nk = options.gamma_nk;
Y = options.Y;
% beta1 = options.beta1;
% KL = options.KL;
w_k = options.w_k;
rho2 = options.rho2;

[N,M] = size(A);
K1 = size(K_all,1);

[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);

N_nk = zeros(N,K1);
for k = 1:K1
    N_nk(:,k) = logmvn(Y, mu_nk(:,:,k), sigma_nk(:,:,:,k), options);
end 

val = -sum(sum(0.5 * gamma_nk .* N_nk));

val = val + calc_prior_sigma(w_k, sigma_jk, K_all, rho2);

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

N_nk = zeros(N,K1);
for k = 1:K1
% parfor k = 1:K1
    N_nk(:,k) = logmvn(Y, mu_nk(:,:,k), sigma_nk(:,:,:,k), options);
end 

val = -sum(sum(0.5 * gamma_nk .* N_nk)) + beta1/2 * trace(A'*KL*A);

%% Calculate derivatives of mu, sigma and A
function der_mu = calc_der_mu(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y, rho1, H, w_k)
[N,M] = size(A);
B = size(mu_jk{1},2);
K1 = size(K_all,1);
K = max(K_all,[],1);

lambda_nk = calc_lambda_nk(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y);
[mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all);

der_mu = mu_jk;
temp = zeros(M,B,K1);
for k = 1:K1
    temp(:,:,k) = A' * lambda_nk(:,:,k) - rho1 * w_k(k) * H * mu_all(:,:,k)';
end

for j = 1:M
    for k = 1:K(j)
        der_mu{j}(k,:) = -sum(temp(j,:,K_all(:,j)==k), 3)';
    end
end


function der_sigma = calc_der_sigma(A, sigma, mu_jk, sigma_jk, K_all, ...
    gamma_nk, Y, rho2, w_k)
[N,M] = size(A);
[~,B] = size(Y);
K = max(K_all,[],1);
K1 = size(K_all,1);
w_jk = w_k2w_jk(w_k,K_all);

[lambda_nk,psi_nk] = calc_lambda_psi_nk(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y);

psi_k = reshape(psi_nk, [B*B,N,K1]);
der_sigma = sigma_jk;
temp = zeros(M,B*B,K1);
for k = 1:K1
    temp(:,:,k) = (A.^2)' * psi_k(:,:,k)';
end

for j = 1:M
    for k = 1:K(j)
        der_sigma{j}(:,:,k) = reshape(-sum(temp(j,:,K_all(:,j)==k), 3)', B, B) ...
              + rho2 * w_jk{j}(k) * sigma_jk{j}(:,:,k);
%             + rho2 * sigma_jk{j}(:,:,k);
    end
end

function der_A = calc_der_A(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, ...
    Y, KL, beta1)
[N,K1] = size(gamma_nk);
M = length(mu_jk);
B = size(mu_jk{1},2);

[lambda_nk,psi_nk] = calc_lambda_psi_nk(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y);

[mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all);
psi_k = reshape(psi_nk, [B*B,N,K1]);

der_A = zeros(N,M);
temp1 = zeros(N,M);
for k = 1:K1
    der_A = der_A - lambda_nk(:,:,k) * mu_all(:,:,k);
    temp1 = temp1 - psi_k(:,:,k)' * sigma_all(:,:,k);
end
der_A = der_A + 2 * A .* temp1;

der_A = der_A + beta1 * KL * A;

%% Calculate lambda_nk and psi_nk
function lambda_nk = calc_lambda_nk(A, sigma, mu_jk, sigma_jk, K_all, gamma_nk, Y)
[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);

[N,B,K1] = size(mu_nk);

lambda_nk = zeros(N,B,K1);
for k = 1:K1
    Y1 = (Y - mu_nk(:,:,k))';
    Y2 = Y1 .* repmat(gamma_nk(:,k)', B, 1);
    sigma_k = block_diag(sigma_nk(:,:,:,k));
    lambda_nk(:,:,k) = reshape(sigma_k \ Y2(:), [B,N])';
end


function [lambda_nk,psi_nk] = calc_lambda_psi_nk(A, sigma, mu_jk, ...
    sigma_jk, K_all, gamma_nk, Y)
[mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, sigma, mu_jk, sigma_jk, K_all);

[N,B,K1] = size(mu_nk);

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

%% Update mu, sigma and A based on the derivative and time step
function mu_jk_new = update_mu(mu_jk, options, delta_t)
der_mu = options.der_mu;

mu_jk_new = mu_jk;

for j = 1:length(mu_jk)
    mu_jk_new{j} = mu_jk{j} - delta_t * der_mu{j};
end

function sigma_jk_new = update_sigma(sigma_jk, options, delta_t)
der_sigma = options.der_sigma;

sigma_jk_new = sigma_jk;

bound_touched = 0;
for j = 1:length(sigma_jk)
    sigma_jk_new{j} = sigma_jk{j} - delta_t * der_sigma{j};
%     sigma_jk_new{j} = sigma_jk{j};
    for k = 1:size(sigma_jk{j},3)
        cov_mat = sigma_jk_new{j}(:,:,k);
        [V,D] = eig((cov_mat + cov_mat')/2);
        d = diag(D);
        if (min(d) < 0)
            bound_touched = bound_touched + 1;
            d(d<0) = 1e-6 * max(abs(d));
            D = diag(d);
        end
        sigma_jk_new{j}(:,:,k) = V*D*V';
    end
end

if bound_touched > 0
%     disp(['There are negative eigenvalues in ', ...
%         num2str(bound_touched), ' matrices of updated sigma_jk']);
end

function A_new = update_A(A, options, delta_t)
der_A = options.der_A;

A_new = A - delta_t * der_A;
A_new = project_to_simplex(A_new);

%% Calculate mu_nk and sigma_nk
function [mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all)
[N,M] = size(A);
K1 = size(K_all,1);
B = size(mu_jk{1},2);

[mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all);

sigma2_I = sigma;

mu_nk = zeros(N,B,K1);
sigma_nk = zeros(B*B,N,K1);
for k = 1:K1
    mu_nk(:,:,k) = A * mu_all(:,:,k)';
    sigma_nk(:,:,k) = sigma_all(:,:,k) * (A.^2)' + repmat(sigma2_I(:),1,N);
end

sigma_nk = reshape(sigma_nk,[B,B,N,K1]);

function [mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all)
mu_all = calc_mu_all(mu_jk, K_all);
sigma_all = calc_sigma_all(sigma_jk, K_all);


function [mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, sigma, mu_jk, sigma_jk, K_all)
[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, sigma, mu_jk, sigma_jk, K_all);

% prec1 = multinv(sigma_nk);
% To avoid out of memory. The running time is similar to the above single
% statement.
prec = zeros(size(sigma_nk));
for k = 1:size(prec,4)
% parfor k = 1:size(prec,4)
    prec(:,:,:,k) = multinv(sigma_nk(:,:,:,k));
end


%% Transform between w_jk and w_k
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
K = K_all2K(K_all);

for j = 1:M
    w_jk{j} = zeros(1,K(j));
end

for j = 1:M
    for l = 1:K(j)
        w = w_k(K_all(:,j)==l);
        w_jk{j}(l) = sum(w);
    end
end

function R = gmm2mean(w_jk, mu_jk)
M = length(w_jk);
B = size(mu_jk{1},2);
R = zeros(M,B);
for j = 1:M
    R(j,:) = sum(repmat(w_jk{j}',1,B) .* mu_jk{j}, 1);
end

function [mu_jk,sigma_jk] = gmm_restore_from_projection(mu_jk,sigma_jk,mapping)
[mu_jk,sigma_jk] = restore_from_projection(mu_jk,sigma_jk,[],mapping.mean,mapping.M);
