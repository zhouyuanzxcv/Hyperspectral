function [ output_args ] = test_scm( input_args )
%TEST_OBJ_FUN Summary of this function goes here
%   Detailed explanation goes here
close all;
addpath('../common');
addpath('../AsterSpectralLibrary');
addpath('../VCA');
addpath('../');

dataset = '00';

options.beta1 = 0;
options.beta2 = 0;
options.rho1 = 0;
options.show_figure = 1;


step_size = 'medium';
fcn_error = @(error_M, error_A) error_M;

disp('------------------------------------------------------------------');
%% generate dataset
switch dataset
    case '00'
        load('toy_image.mat');
        
%         sigmas = noise_est_roger(I);
%         disp('Roger error')
%         mdiff(Y_noise,sigmas)
                
        M = 6;
        
        options.beta1 = 10;
        options.beta2 = 10;
        options.rho1 = 0.001; % Set it to 0.1 to see uncertainty
        
        options.skip_uncertainty = 0;
        options.use_single_variance_for_unmixing = 0;
        options.use_single_variance_for_uncertainty = 0;

        step_size = [10, 10, 10];
        fcn_error = @(error_M, error_A) error_M*10 + error_A;
    case '01'
        load('toy_image_snr_40.mat');
        R_gt = E; A_gt = A;
                
        M = 4;
        
        options.beta1 = 10;
        options.beta2 = 10;
        options.rho1 = 0.1; % Shows uncertainty
        
        options.use_single_variance_for_uncertainty = 0;
    case '2'
        load('../../data/PaviaUniversity_corrected.mat');
        % remove gravel and bitumen
        R_gt([3,7],:) = []; 
        names([3,7]) = [];
        A_gt(:,:,[3,7]) = [];
        A_gt = double(A_gt);
                
        M = 7;
%         options.beta1 = 0.01;
%         options.beta2 = 0.02;
%         options.rho1 = 0.05;   

        options.beta1 = 5;
        options.beta2 = 10;
        options.rho1 = 0.05;
        
   case '61'
        load('muufl_gulfport_roi_corrected.mat');
        figure, imshow(rgb);
        figure, imshow(I_gt,[1,1,1;gt_colors]);
        colorbar('YTickLabel',tick_labels,'YDir','reverse','YTick',...
            (0.5:1:size(A_gt,3)+0.5)','YLim',[0 size(A_gt,3)]);

        A_gt = double(A_gt);
        
        M = 6;

        options.beta1 = 10;
        options.beta2 = 50;
        options.rho1 = 0.01;
    otherwise
        disp('Unknown method to generate dataset.')
end


if isempty(names)
    names = cell(1,M);
    for i = 1:M
        names{i} = ['endmember ',num2str(i)];
    end
end

%% Run SCM
[Y,A_gt,rows,cols] = reshape_hsi(I,A_gt);
[N,B] = size(Y);

if exist('ws_gt','var') && ~isempty(ws_gt)
    endmember_scatter_plot_end_var(Y,ws_gt,mus_gt,sigmas_gt,names);
    R_gt = zeros(M,B);
    for j = 1:M
        R_gt(j,:) = mean(Y(A_gt(:,j)==1,:), 1);
    end
else
    endmember_scatter_plot(Y,R_gt,names);
end
set(gcf,'name','Scatter plot of the original image with ground truth');

if ~isempty(A_gt)
    SNR = calc_snr(Y,A_gt,R_gt);
    disp(['SNR: ',num2str(SNR), 'dB']);
end

if 0 % Set to 1 to automatically determine the best parameters
    options.show_figure = 0;
    options.skip_uncertainty = 1;

    pa = ParameterAnalysis();
    list_params = {'beta1','beta2','rho1'};
    range = [1e-6,1e6;1e-6,1e6;1e-6,1e6]';
    fcn_run = @(options) try_algo_for_param(I,M,options,R_gt,A_gt, ...
        names,wl,fcn_error);
    options = pa.autoParamSelection(fcn_run, options, list_params, ...
        step_size, range);
    
    options.skip_uncertainty = 0;
    options.show_figure = 1;
end

[A,R,mu,sigma,var_dirs,var_amts] = scm(I,M,options);

options1 = [];
options1.var_dirs = var_dirs;
options1.var_amts = var_amts;

if M <= 6, plot_row = ceil(M/2); else plot_row = 2; end
options1.plot_row = plot_row;

[error_M,error_A,best_p] = compare_2_endmembers(R_gt,R,A_gt,A, ...
    rows,cols,names,wl,1,options1);

if ndims(sigma) == 3
    % permute the sigmas (the covariance matrices)
    sigma = sigma(:,:,best_p*(1:M)');

    var_dirs = best_p * var_dirs;
    var_amts = best_p * var_amts;
    
    % plot the variations
    figure,plot(var_dirs');
    disp(['uncertainty amounts: ', num2str(var_amts')]);
    disp(['total uncertainty: ', num2str(sum(var_amts))]);
elseif length(sigma) > 1
    sigma = best_p*sigma;
    disp(['sigma: ', num2str(sigma')]);
end

A = A*best_p';
R = best_p*R;

save('result_scm.mat','A','R','mu','sigma','var_dirs','var_amts');

% show pixel cloud
endmember_scatter_plot(Y,R,names);


function val = try_algo_for_param(I, M, options, R_gt, A_gt, names, wl, fcn_error)
close all;
[rows,cols,B] = size(I);
[A,R,mu,sigma,var_dirs,var_amts] = scm(I,M,options);
[error_M,error_A,best_p] = compare_2_endmembers(R_gt, R, A_gt, A, ...
    rows,cols,names,wl,0);

val = fcn_error(error_M, error_A);