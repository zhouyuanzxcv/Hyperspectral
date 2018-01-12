function [A,R,w_jk,mu_jk,sigma_jk,extra] = gmm_hu_ex(I, endmembers, options)
%GMM_EX Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    options = [];
end

if max(I(:)) > 10 || max(endmembers{1}(:)) > 10
    disp('Warning! I or endmember spectra are not in the range 0 - 1.');
end

if class(I) ~= class(endmembers{1})
    disp('Error! I and endmember spectra do not have the same data type');
end

options = insert_param_when_absent(options, 'reduced_dim', 10);
options = insert_param_when_absent(options, 'max_num_comp', 4);
options = insert_param_when_absent(options, 'show_fig', 0);

options = insert_param_when_absent(options, 'fix_w_k', 1);
options = insert_param_when_absent(options, 'fix_mu_jk', 1);
options = insert_param_when_absent(options, 'fix_sigma_jk', 1);

use_predefined_projection_components = parse_param(options, ...
    'use_predefined_projection_components', 0);

[Y_ori,~,rows,cols] = reshape_hsi(I,[]);
M = length(endmembers);

if use_predefined_projection_components
    mapping = options.project_mapping;
    Y = gmm_project(Y_ori, mapping);
    components = options.predefined_components;
    K = components.K;
    w_jk = components.w_jk;
    mu_jk = components.mu_jk;
    sigma_jk = components.sigma_jk;
else
    % project_mode = 'image' or 'endmembers'
    project_mode = parse_param(options,'project_mode','endmembers');
    
    if strcmp(project_mode, 'image')
        [Y, mapping] = pca(Y_ori, options.reduced_dim);
    elseif strcmp(project_mode, 'endmembers')
        mapping = calc_projection_from_library(endmembers, options.reduced_dim);        
        Y = gmm_project(Y_ori, mapping);
    elseif strcmp(project_mode, 'combined')
        % Gives pretty bad results, should figure out a new implemenation
        mapping = calc_projection_from_image_library(Y_ori, endmembers, options.reduced_dim);
        Y = gmm_project(Y_ori, mapping);
    elseif strcmp(project_mode, 'custom')
        mapping = options.project_mapping;
        Y = gmm_project(Y_ori, mapping);
    else
        disp('ERROR: Unknown projection mode');
    end
    options = insert_param_when_absent(options, 'project_mapping', mapping);
    
    % Y = gmm_project(Y_ori, mapping);
    [K,w_jk,mu_jk,sigma_jk] = estimate_components(endmembers, mapping, options);
end

if options.show_fig
    names = options.names;
    endmember_scatter_plot_end_var(Y,w_jk,mu_jk,sigma_jk,names);
end

options.w_jk = w_jk;
options.mu_jk = mu_jk;
options.sigma_jk = sigma_jk;
options.K = K;

options = insert_param_when_absent(options, 'beta1', 0);
options = insert_param_when_absent(options, 'beta2', 0);
% options = insert_param_when_absent(options, 'convergence_thresh', 0.002);

% no need to project back to the original dimension, gmm_hu will do it.
[A,R,w_jk,mu_jk,sigma_jk,extra] = gmm_hu(I,M,options);

function mapping = calc_projection_from_image_library(Y_ori, endmembers, reduced_dim)
spectra_train = cell2mat(endmembers');
if size(spectra_train,1) > size(Y_ori,1)
    Y_ori = repmat(Y_ori, round(size(spectra_train,1) / size(Y_ori,1)), 1);
end
spectra_train = cat(1, spectra_train, Y_ori);
[~,mapping] = pca(spectra_train, reduced_dim);


function [beta1,beta2] = automatic_params(I,M,options)
init_beta1 = 10;
init_beta2 = 10;

options.fix_w_k = 0;
options.convergence_thresh = 0.1;

w_jk = options.w_jk;

beta1s = 10.^(-5:5)';
beta2s = 10.^(-5:5)';

errs = Inf(length(beta1s), length(beta2s));
% A_errs = Inf(length(beta1s), length(beta2s));

[~,ind1] = min(abs(beta1s - init_beta1));
[~,ind2] = min(abs(beta2s - init_beta2));
while 1
    ind1s = [ind1-1, ind1, ind1+1];
    ind2s = [ind2-1, ind2, ind2+1];
%     ind2s = ind2;
    for i = ind1s
        for j = ind2s
            % only the inner errors will be updated so the current point
            % will not move to the boundary
            if i > 1 && i < length(beta1s) && j > 1 && ...
                    j < length(beta2s) && errs(i,j) == Inf
                options.beta1 = beta1s(i);
                options.beta2 = beta2s(j);
                [A,R,w_jk1,mu_jk,sigma_jk] = gmm_hu(I,M,options);
                title = ['beta1: ',num2str(options.beta1),', beta2: ',num2str(options.beta2)];
                show_abundances(A, rows, cols, title);
%                 A_errs(i,j) = calc_abundance_error(A_gt,A);
                errs(i,j) = mdiff(w_jk,w_jk1);
            end
        end
    end
    if errs(ind1,ind2) <= errs(ind1-1,ind2) && errs(ind1,ind2) <= errs(ind1+1,ind2) ...
            && errs(ind1,ind2) <= errs(ind1,ind2-1) && errs(ind1,ind2) <= errs(ind1,ind2+1)
        break;
    else
        errs1 = errs(ind1-1:ind1+1,ind2-1:ind2+1);
        [~,ind] = min(errs1(:));
        [y,x] = ind2sub([3 3], ind);
        ind1 = ind1 + y - 2;
        ind2 = ind2 + x - 2;
    end
end

beta1 = beta1s(ind1);
beta2 = beta2s(ind2);


