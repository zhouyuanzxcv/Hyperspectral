function [ output_args ] = test_gmm_ex( input_args )
%TEST_GMM_HU_EX Summary of this function goes here
%   Detailed explanation goes here
close all;

dataset = '7';
[endmembers,I,Y,R_gt,A_gt,names,wl] = prepare_supervised_unmixing(dataset);
[rows,cols,B] = size(I);

D = 0.001^2 * eye(B);
options.beta1 = 0;
options.beta2 = 0;
options.show_fig = 1;
options.names = names;
options.D = D;
options.project_mode = 'image';
options.convergence_thresh = 0.0001;

[A,R,w_jk,mu_jk,sigma_jk,extra] = gmm_hu_ex(I, endmembers, options);
E = gmm_hu_endmember(I,A,D,w_jk,mu_jk,sigma_jk);
save('result_gmm.mat','A','R','E','w_jk','mu_jk','sigma_jk');

mdiff(A,A_gt);
show_abundances(A,rows,cols);
replay_scatter_abund(extra.frames_scatter, extra.frames_abund);

end

