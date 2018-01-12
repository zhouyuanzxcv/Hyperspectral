function [ output_args ] = test_ncm(dataset, use_pca)
%TEST_GMM Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
    dataset = '001';
end
if nargin < 2
    use_pca = 1;
end

close all

%% load dataset
[endmembers,I,Y,R_gt,A_gt,names,wl] = prepare_supervised_unmixing(dataset);

I1 = retrieve_rgb(I,wl);
figure('name','RGB image of the original image');
imshow(I1);
[rows,cols,B] = size(I);
[N,B] = size(Y);
M = size(A_gt,2);

D = 0.001^2 * eye(B);

options.convergence_thresh = 0.0001;
options.D = D;
options.A_gt = A_gt; % for testing

options.max_num_comp = 1;
options.beta1 = 0;
options.beta2 = 0;
options.names = names;

if use_pca
    options.project_mode = 'image';
else
    mapping = struct('mean',zeros(1,B),'M',eye(B));
    options.project_mode = 'custom';
    options.project_mapping = mapping;
end

[A,R,w_jk,mu_jk,sigma_jk] = gmm_hu_ex(I, endmembers, options);

E = gmm_hu_endmember(I,A,D,w_jk,mu_jk,sigma_jk);

% mdiff(A,A_gt);

%% Permute the results to accord with the GT for comparison
[error_M,error_A,best_p] = compare_2_endmembers(R_gt, R, A_gt, A, ...
    rows,cols,names,wl,1);

A = A*best_p';
R = best_p*R;
for i = 1:size(E,3)
    E(:,:,i) = best_p * E(:,:,i);
end
w_jk = w_jk(best_p*(1:M)');
mu_jk = mu_jk(best_p*(1:M)');
sigma_jk = sigma_jk(best_p*(1:M)');

save(['result_ncm_',num2str(use_pca),'_',dataset,'.mat'],'A','R','E',...
    'w_jk','mu_jk','sigma_jk');

opts = struct('show_approx',1,'w_jk',{w_jk},'mu_jk',{mu_jk},'sigma_jk',...
    {sigma_jk},'legend_names',{{'Estimated distribution'}});
hist_end_var(Y,A_gt,names,1,opts);

