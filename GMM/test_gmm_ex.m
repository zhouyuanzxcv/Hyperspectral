function [ output_args ] = test_gmm_ex( input_args )
%TEST_GMM_HU_EX Summary of this function goes here
%   Detailed explanation goes here
close all;

dataset = '7';
[endmembers,I,Y,R_gt,A_gt,names,wl] = prepare_supervised_unmixing(dataset);
[rows,cols,B] = size(I);

% noise covariance matrix. it can be zero.
D = 0.001^2 * eye(B); 

% smoothness and sparsity constraints on the abundances
options.beta1 = 0;
options.beta2 = 0;

% show intermediate results (1) or not (0)
options.show_fig = 1;

options.names = names;
options.D = D;

% project_mode can be
%   'image' - apply PCA on the pixels of the image
%   'endmembers' - concatenate the spectra in the library and use PCA on them
options.project_mode = 'image';

% threshold of convergence
options.convergence_thresh = 0.0001;

% calculate abundance maps
[A,R,w_jk,mu_jk,sigma_jk,extra] = gmm_hu_ex(I, endmembers, options);

% calculate endmembers per pixel
E = gmm_hu_endmember(I,A,D,w_jk,mu_jk,sigma_jk);

save('result_gmm.mat','A','R','E','w_jk','mu_jk','sigma_jk');

mdiff(A,A_gt);
show_abundances(A,rows,cols);
replay_scatter_abund(extra.frames_scatter, extra.frames_abund);

end

