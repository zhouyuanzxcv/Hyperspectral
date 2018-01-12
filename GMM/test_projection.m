function [ output_args ] = test_projection( input_args )
%TEST_PROJECTION Summary of this function goes here
%   Detailed explanation goes here
%% show scatter plots
load('toy_image_end_var_snr60_1231.mat');
[Y_gt,A_gt,rows,cols] = reshape_hsi(I,A_gt);
[N,B] = size(Y_gt);
M = size(A_gt,2);
endmember_scatter_plot_end_var(Y_gt,ws_gt,mus_gt,sigmas_gt,names);
set(gcf,'name','Scatter plot of the original image with ground truth');

reduced_dim = 10;
[Y, mapping] = pca(Y_gt, reduced_dim);
[mus_gt1,sigmas_gt1] = project2ortho(mus_gt,sigmas_gt,mapping.mean,mapping.M);
endmember_scatter_plot_end_var(Y,ws_gt,mus_gt1,sigmas_gt1,names);
R_gt = zeros(M,B);
R = zeros(M,reduced_dim);
for j = 1:M
    R_gt(j,:) = mean(Y_gt(A_gt(:,j)==1,:), 1);
    R(j,:) = mean(Y(A_gt(:,j)==1,:), 1);
end
set(gcf,'name','Scatter plot of the dimensionality-reduced image');

%% show histograms
opts = struct('show_approx',1,'w_jk',{ws_gt},'mu_jk',{mus_gt},'sigma_jk',...
    {sigmas_gt},'legend_names',{{'Ground truth GMM'}});
hist_end_var(Y_gt,A_gt,names,1,opts);
set(gcf,'name','Histogram of the ground truth pure pixels in the original dimension');

opts.mu_jk = mus_gt1;
opts.sigma_jk = sigmas_gt1;
hist_end_var(Y,A_gt,names,1,opts);
set(gcf,'name','Histogram of the ground truth pure pixels in the reduced dimension');

[mu_jk1,sigma_jk1,R1] = restore_from_projection(mus_gt1,sigmas_gt1,R,mapping.mean,mapping.M);

% test positive definiteness of covariance matrices
for j = 1:M
    for k = 1:size(sigma_jk1,3)
        [R,err] = cholcov(sigma_jk1{j}(:,:,k),0);
        if err ~= 0
            disp(['The ',num2str(k),'th component of the ',num2str(j),'th endmember is not SPD']);
        end
    end
end
% test difference from the original means and covariances
disp('Mu reconstruction error: '); mdiff(mus_gt,mu_jk1)
disp('Sigma reconstruction error: '); mdiff(sigmas_gt,sigma_jk1)

opts.mu_jk = mu_jk1;
opts.sigma_jk = sigma_jk1;
hist_end_var(Y_gt,A_gt,names,1,opts);
