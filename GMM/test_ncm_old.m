function [ output_args ] = test_ncm( input_args )
%TEST_GMM Summary of this function goes here
%   Detailed explanation goes here

close all

%% load dataset

ws_gt = [];

options = struct('reduced_dim',10,'max_num_comp',1);

dataset = '21';
switch dataset
    case '02'
        load('toy_image_end_var_snr60_1231.mat');
        M = 4;
        options.eta = 0.05;
        options.Y_noise = 1e-4;

        options.beta1 = 0.2;
        options.beta2 = 0.01;
        options.rho1 = 0;
        options.rho2 = 0;
        options.sigma0 = 1e-2;
        options.beta2_decay = 1;
    case '21'
        load('../../Data/PaviaUniversity_A.mat');
        A_gt = double(A_gt);
        M = 5;
%         reduced_dim = 10;
        
        options.eta = 0.05;
        options.Y_noise = 1e-4;

        options.beta1 = 5;
        options.beta2 = 5; % 2e4
        options.rho1 = 0;
        options.rho2 = 0;
        options.sigma0 = 0.1;
        options.shrink_size = 2;
        options.beta2_decay = 0.05;
    case '62'
        load('../../data/muufl_gulfport_B.mat');
        A_gt = double(A_gt);
        M = 5;
        
        options.eta = 0.05;
        options.Y_noise = 1e-3;

        options.beta1 = 5;
        options.beta2 = 5; % 2e4
        options.rho1 = 0;
        options.rho2 = 0;
        options.sigma0 = 0.1;
        options.shrink_size = 1;
        options.beta2_decay = 0.02;

%         options.convergence_thresh = 0.001;
    otherwise
end
I1 = retrieve_rgb(I,wl);
figure('name','RGB image of the original image');
imshow(I1);

[Y,A_gt,rows,cols] = reshape_hsi(I,A_gt);
[N,B] = size(Y);

if ~isempty(ws_gt)
    endmember_scatter_plot_end_var(Y,ws_gt,mus_gt,sigmas_gt,names);
    R_gt = zeros(M,B);
    for j = 1:M
        R_gt(j,:) = mean(Y(A_gt(:,j)==1,:), 1);
    end
else
    endmember_scatter_plot(Y,R_gt,{'Ground Truth'});
end
set(gcf,'name','Scatter plot of the original image with ground truth');

if ~isempty(ws_gt)
    opts = struct('show_approx',1,'w_jk',{ws_gt},'mu_jk',{mus_gt},'sigma_jk',...
        {sigmas_gt},'legend_names',{{'Ground truth GMM'}});
else
    opts = struct('legend_names',{{}});
end
hist_end_var(Y,A_gt,names,1,opts);
set(gcf,'name','Histogram of the ground truth pure pixels vs GMM');

I1 = reshape(Y, [rows, cols, B]);

[A,R,w_jk,mu_jk,sigma_jk] = gmm_hu(I1,M,options);

%% Permute the results to accord with the GT for comparison
[error_M,error_A,best_p] = compare_2_endmembers(R_gt, R, A_gt, A, ...
    rows,cols,names,wl,1);

A = A*best_p';
R = best_p*R;
w_jk = w_jk(best_p*(1:M)');
mu_jk = mu_jk(best_p*(1:M)');
sigma_jk = sigma_jk(best_p*(1:M)');

save('result_ncm.mat','A','R','w_jk','mu_jk','sigma_jk');

opts = struct('show_approx',1,'w_jk',{w_jk},'mu_jk',{mu_jk},'sigma_jk',...
    {sigma_jk},'legend_names',{{'Estimated distribution'}});
hist_end_var(Y,A_gt,names,1,opts);

