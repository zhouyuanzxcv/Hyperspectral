function [A,R,mu,sigma,var_dirs,var_amts] = scm(I,M,options)
%SCM Spatial compositional model
% Input:
%   I: row*col*B image data
%   M: number of endmembers
%   options: structure for additional parameters
% Output:
%   A: abundances (N by M)
%   R: endmembers (M by B)
%   mu: noise std
%   sigma: uncertainty of endmembers
%   var_dirs: uncertainty directions
%   var_amts: uncertainty amounts
%
% Note: options is a structure containing eta,beta1,beta2,rho1,...
%   show_figure,init_mode, etc.
% 
start = tic;

if ndims(I) == 3
    [rows, cols, B] = size(I);
    N = rows*cols;
    Y = reshape(I, rows*cols, B);
else
    Y = I;
    [N,B] = size(I);
    I = reshape(I,[N 1 B]);
end

% set parameters
beta1 = 10;
beta2 = 10;
rho1 = 0.01;
rho2 = 0;
eta = 0.05;
show_figure = 1;
skip_uncertainty = 0;
init_mode = 2;
seg_map = [];

iter_num = 300;

use_single_variance_for_unmixing = 0;
use_single_variance_for_uncertainty = 0;

use_fix_R = 0;

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


% beta1 = beta1*B;
% beta2 = beta2*B;
beta1 = beta1*B/M;
beta2 = beta2*B/M;
rho1 = rho1*N/M^2;
rho2 = rho2*N/M;

if beta1 < 1e-9
    beta1 = 1e-9;
end

[K,H] = calc_Laplacians(I,seg_map,M,eta,beta1,beta2);

[A,R] = scm_init(Y, M, init_mode, options);
if show_figure
    endmember_scatter_plot(Y,R,names);
    set(gcf,'name','Scatter plot of the initialization');
end

if use_fix_R
    R = options.fix_R;
end

disp(['Initial endmember closeness: ',num2str(trace(R'*H*R))]);

options.rows = rows;
options.cols = cols;

if use_single_variance_for_unmixing
    % When the variances of noise of each band are equal
    [A,R] = iterate_single_variance(Y,A,R,K,H,beta1,rho1,...
        iter_num,names,show_figure,options);
else
    [A,R] = iterate_distinct_variances(Y,A,R,K,H,beta1,rho1,...
        iter_num,names,show_figure,options);
end

if ~skip_uncertainty
    if use_single_variance_for_uncertainty
        % Estimate the uncertainty
        mu = sqrt(sum(sum((Y-A*R).^2))/(N*B));
        sigma0 = 0.1;
        sigma_max = 1;
        [sigma,mu,var_dirs,var_amts] = scm_ex(Y,A,R,sigma0,mu,sigma_max);
    else
        % When the variances of noise are distinct
        [d,S] = solve_for_d_S(Y,A,R,options);
        [sigma,var_dirs,var_amts] = calc_uncertainty_range(S,d);
        mu = 1./d;
    end
else
    sigma = [];
    var_dirs = [];
    var_amts = [];
    mu = 0;
end

disp(['Elapsed time for SCM: ',num2str(toc(start))]);


function [A,R] = iterate_single_variance(Y,A,R,K,H,beta1,rho1,...
    iter_num,names,show_figure,options)
delta_t0 = 1e-2;
rows = options.rows;
cols = options.cols;
B = size(Y,2);

d = ones(B,1);

errors = calc_error(Y,A,R,d,K,H,beta1,rho1);

epc = EvolveProcessCapture();
epc = epc.begin(iter_num);

for iter = 1:iter_num
    % solve for A
    obj_fun_A = @(A) sum(sum((Y-A*R).^2)) + beta1*trace(A'*K*A);
    fun_der = @(A) -Y*R' + A*R*R' + beta1*K*A;

    A = projected_gradient_descent(A, obj_fun_A, fun_der, ...
        @project_to_simplex, delta_t0);
    % solve for R
    if isfield(options,'use_fix_R') && options.use_fix_R
        R = options.fix_R;
    else
        R = inv(A'*A + rho1*H) * A' * Y;
    end

    % calc error
    err = calc_error(Y,A,R,d,K,H,beta1,rho1);
    errors = [errors;err];
    
    if show_figure
        epc = epc.updateEndmFixed(Y,A,rows,cols,R,names,iter);
    end

    if test_convergence(errors, 1e-5), break; end
    
    if mod(iter,5) == 0
        disp(['Process iteration ',num2str(iter)]);
    end
end

epc = epc.end(iter);

if show_figure
    figure, plot([1:length(errors)]',errors);
    xlabel('Iteration number');
    ylabel('Energy');
end


function [A,R] = iterate_distinct_variances(Y,A,R,K,H,beta1,rho1,...
    iter_num,names,show_figure,options)
% iter_num = 500;
rows = options.rows;
cols = options.cols;
[N,B] = size(Y);
M = size(A,2);

d = sqrt(mean((Y-A*R).^2,1)');
d = 1./d;
% d = ones(B,1);

epc = EvolveProcessCapture();
epc = epc.begin(iter_num);

disp('Start estimating A and R');

errors = calc_error(Y,A,R,d,K,H,beta1,rho1);
for iter = 1:iter_num
    % solve for A
    obj_fun_A = @(A) eval_obj_fun_A(Y,A,R,d,K,beta1);
    fun_der_A = @(A) calc_der_A(Y,A,R,d,K,beta1);

    A = projected_gradient_descent(A, obj_fun_A, fun_der_A, ...
        @project_to_simplex, 1e-3);
    
    % solve for R
    invA1A = inv(A'*A + rho1*H);
    R = invA1A*A'*Y;
%     R = solve_for_R(Y,A,H,d,rho1);

    % solve for D
%     d = sqrt(sum((Y-A*R).^2,1)'/N);
    d = sqrt((sum((Y-A*R).^2,1)' + rho1*diag(R'*H*R))/N);
    d = 1./d;
        
    % calc error
    err = calc_error(Y,A,R,d,K,H,beta1,rho1);
    errors = [errors;err];
    if test_convergence(errors, 1e-5), break; end
    
    if mod(iter,10) == 0
        if show_figure
            epc = epc.updateEndmFixed(Y,A,rows,cols,R,names,iter);
        end
        disp(['Process iteration ',num2str(iter),'. The objective ',...
            'function has value ',num2str(err)]);
    end
end

epc = epc.end(iter);


if show_figure
    figure, plot([1:length(errors)]',errors);
    xlabel('Iteration number');
    ylabel('Energy');
end


function der = calc_der_A(Y,A,R,d,K,beta1)
M = size(R,1);
R1 = repmat(d.^2,[1,M]) .* (R');
der = -Y*R1 + A*R*R1 + beta1*K*A;

der_norm = max(abs(der(:)));
der = der/der_norm;


function val = eval_obj_fun_A(Y,A,R,d,K,beta1)
N = size(Y,1);
YAR1 = repmat(d,[1,N]) .* (Y-A*R)';
val = sum(sum(YAR1.^2)) + beta1*trace(A'*K*A);


function R = solve_for_R(Y,A,H,d,rho1)
M = size(A,2);
B = length(d);
R = zeros(M,B);
C = A'*A;
F = A'*Y;
for k = 1:B
    R(:,k) = (C + rho1*(d(k)^(-2))*H) \ F(:,k);
end


function err = calc_error(Y,A,R,d,K,H,beta1,rho1)
[N,B] = size(Y);
YAR = Y-A*R;
lsq_err = sum(sum((repmat(d,[1,N]) .* YAR').^2));
log_D = -N * sum(log(d.^2));
A_err = beta1 * trace(A'*K*A);
R_close = rho1 * trace(repmat(d.^2,[1,B]) .* (R'*H*R));
% R_close = rho1 * trace(R'*H*R);

err = lsq_err + log_D + A_err + R_close;


