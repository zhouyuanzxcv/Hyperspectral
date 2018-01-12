function [ output_args ] = test_gmm(dataset)
%TEST_GMM Summary of this function goes here
%   Detailed explanation goes here
addpath('../NCM');

close all

% select dataset - '21' for Pavia University, '62' for the gulfport dataset
% '00' for generating unsupervised synthetic dataset
% '001' or '002' for running on the unsupervised synthetic image 1 or 2
if nargin < 1
    dataset = '001';
end

%% load dataset
SNR = 60;
% asphalt, basalt, concrete, conifer, grass, limestone, quartzite
chosen = [6,2,3,1];
rows = 60;
cols = 60;

options = struct('reduced_dim',10);
options.convergence_thresh = 0.0001;
options.show_fig = 0;

switch dataset
    case '00'
        Y_noise = 1e-3;
        Y_noise = rand(200,1) * Y_noise;

        gen_im_opts = [];
        gen_im_opts.chosen = chosen;
        gen_im_opts.noise_sigma = Y_noise;
        
        if 0 % synthetic unsupervised A
%             [I,E_gt,A_gt,names,wl,ws_gt,mus_gt,sigmas_gt,Y_noise] = ...
%                 generate_toy_image_end_var(rows, cols, Y_noise, chosen, 1);
            [I,E_gt,A_gt,names,wl,Y_noise,extra] = generate_toy_image_unified(...
                rows, cols, 'gmm', 'quadrant', 0, gen_im_opts);
            ws_gt = extra.ws;
            mus_gt = extra.mus;
            sigmas_gt = extra.sigmas;
            save('toy_image_end_var_3.mat','I','E_gt','A_gt','names','wl',...
                'ws_gt','mus_gt','sigmas_gt','Y_noise');
        elseif 0 % synthetic unsupervised B
            gen_im_opts.material_gaussian_centers_normalized = [0.5,0.25;...
                0.25,0.75; 0.75,0.75]; % (x,y)
            gen_im_opts.material_gaussian_sigma_normalized = 0.25;
            gen_im_opts.material_num_samples = 150;
            gen_im_opts.material_shape_width = 1.5;
            gen_im_opts.material_shape_width_sd = 0.5;
            gen_im_opts.material_gaussian_coefficient = 10;

            [I,E_gt,A_gt,names,wl,Y_noise,extra] = generate_toy_image_unified(...
                rows, cols, 'gmm', 'gaussian_dotted', 0, gen_im_opts);
            ws_gt = extra.ws;
            mus_gt = extra.mus;
            sigmas_gt = extra.sigmas;
            save('toy_image_end_var_3_B.mat','I','E_gt','A_gt','names','wl',...
                'ws_gt','mus_gt','sigmas_gt','Y_noise');
        else
%             load('toy_image_end_var_3.mat');
            load('toy_image_end_var_3_B.mat');
        end
        
%         noise_sigmas = noise_est_roger(I);
%         disp('The error of noise estimation is ');
%         mdiff(Y_noise,noise_sigmas);
%         
%         noise_sigmas = noise_est_mlr(I);
%         mdiff(Y_noise,noise_sigmas);
        
        M = 4;
        options.beta1 = 5;
        options.beta2 = 5; % 0.005
        options.shrink_size = 5;
        options.beta2_decay = 0.05;
    case '001'
        load('toy_image_end_var_3.mat');

        M = 4;
        options.beta1 = 5;
        options.beta2 = 5; % 0.005
        options.shrink_size = 5;
        options.beta2_decay = 0.05;
    case '002'
        load('toy_image_end_var_3_B.mat');

        M = 4;
        options.beta1 = 0.1;
        options.beta2 = 0.1; % 0.005
        options.shrink_size = 5;
        options.beta2_decay = 0.05;
    case '01'
        load('toy_image_end_var_snr60_1111.mat');
        M = 4;
        options.beta1 = 2000;
        options.beta2 = 50; % 2e4
        options.eta = 0.05;
        options.rho1 = 5;
    case '2'
%         [I,R_gt,A_gt,names,wl] = load_pavia_university;
        load('../../data/PaviaUniversity_corrected.mat');
        % remove gravel and bitumen
        R_gt([3,7],:) = []; 
        names([3,7]) = [];
        if 1
            A_gt(:,:,[3,7]) = [];
        else
            to_be_merged = [8,1];
            merged = [3,7];
            for j = 1:2
                A1 = A_gt(:,:,to_be_merged(j));
                A2 = A_gt(:,:,merged(j));
                A3 = A2 + A1;
                A_gt(:,:,to_be_merged(j)) = A3;
            end
            A_gt(:,:,[3,7]) = [];
        end
        M = 7;

        options.beta1 = 0.5;
        options.beta2 = 0.5; % 2e4
        options.rho1 = 0.01;
        options.rho2 = 0;
        options.shrink_size = 2;
        options.beta2_decay = 0.5;

    case '21'
        load('../../data/PaviaUniversity_A.mat');
        A_gt = double(A_gt);
        M = 5;
        
        options.beta1 = 5;
        options.beta2 = 5; % 2e4
        options.rho1 = 0;
        options.rho2 = 0;
        options.shrink_size = 2;
        options.beta2_decay = 0.05;        
    case '3'
        [I,R_gt,A_gt,names,wl] = load_neon;
    case '6'
        [I,R_gt,A_gt,names,wl] = load_gulfport;
    case '62'
        load('../../data/muufl_gulfport_B.mat');
        M = 5;
        
%         remove_bands = [1:5,72-9:72];
%         [I,R_gt,wl] = remove_noisy_bands(I,R_gt,wl,remove_bands);
        
        options.beta1 = 5;
        options.beta2 = 5; % 2e4
        options.rho1 = 0;
        options.rho2 = 0;
        options.shrink_size = 1;
        options.beta2_decay = 0.05;
    otherwise
end
if exist('rgb','var') && ~isempty(rgb)
    I1 = rgb;
else
    I1 = retrieve_rgb(I,wl);
end
figure('name','RGB image of the original image');
imshow(I1);

[Y,A_gt,rows,cols] = reshape_hsi(I,A_gt);
[N,B] = size(Y);

if exist('ws_gt','var') && ~isempty(ws_gt)
    endmember_scatter_plot_end_var(Y,ws_gt,mus_gt,sigmas_gt,names);
    R_gt = zeros(M,B);
    for j = 1:M
        R_gt(j,:) = mean(Y(A_gt(:,j)==1,:), 1);
    end
else
    endmember_scatter_plot(Y,R_gt,{'Ground Truth'});
end
set(gcf,'name','Scatter plot of the original image with ground truth');

if exist('ws_gt','var') && ~isempty(ws_gt)
    opts = struct('show_approx',1,'w_jk',{ws_gt},'mu_jk',{mus_gt},'sigma_jk',...
        {sigmas_gt},'legend_names',{{'Ground truth GMM'}});
else
    opts = struct('legend_names',{{}});
end
hist_end_var(Y,A_gt,names,1,opts);
set(gcf,'name','Histogram of the ground truth pure pixels vs GMM');

I1 = reshape(Y, [rows, cols, B]);

D = 0.001^2 * eye(B);
options.D = D;

% [N,B] = size(Y);
% options.A_gt = A_gt;
% options.fix_A = 1;
% options.fix_sigma_jk = 1;

[A,R,w_jk,mu_jk,sigma_jk,extra] = gmm_hu(I1,M,options);
E = gmm_hu_endmember(I1,A,D,w_jk,mu_jk,sigma_jk);


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

if options.show_fig % Play the movie 
    replay_scatter_abund(extra.frames_scatter, extra.frames_abund);
end

save(['result_gmm_',dataset,'.mat'],'A','R','E','w_jk','mu_jk','sigma_jk');

opts = struct('show_approx',1,'w_jk',{w_jk},'mu_jk',{mu_jk},'sigma_jk',...
    {sigma_jk},'legend_names',{{'Estimated distribution'}});
hist_end_var(Y,A_gt,names,1,opts);


function [ws,mus,sigmas] = create_endmember_params_1_1_1_1(E)
[M,B] = size(E);
ws = cell(1,M);
mus = cell(1,M);
sigmas = cell(1,M);

for j = 1:M
    ws{j} = 1;
    mus{j} = E(j,:);
    sigmas{j} = 0.01^2*eye(B);
end

% endmember variability of basalt
idx = 2;
major_dir = ones(B,1);
major_dir = major_dir/sqrt(sum(major_dir.^2));
sigmas{idx} = 0.01^2*eye(B) + 0.05^2*major_dir*major_dir';

% endmember variability of concrete
idx = 3;
major_dir = ones(B,1);
major_dir = major_dir/sqrt(sum(major_dir.^2));
sigmas{idx} = 0.01^2*eye(B) + 0.03^2*major_dir*major_dir';