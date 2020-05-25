function [error,total_time,As] = run_end_var_algos(algo, dataset, options)
%RUN_END_VAR_ALGOS Summary of this function goes here
%   Detailed explanation goes here
addpath('../common');
addpath('../GMM');
addpath('competing_methods/BCM');
addpath('competing_methods/AAM');
addpath('competing_methods/MESMA');
addpath('competing_methods/unmixP_NCM');

addpath('process_envi_data');
addpath('analysis');

if nargin < 1
    algo = 'GMM';
end
if nargin < 2
    % the format is 'sba_im_#_lib_#' where '#' is either 4 or 16,
    % indicating the spatial resolution for the image or library
    dataset = 'sba_im_16_lib_16';
end
if nargin < 3
    options = [];
%     options.select_image_indices = [1,2,9,13,15,17,27,33,35,39,58,62,64];
end

select_image_indices = parse_param(options, 'select_image_indices', []);
savefile_postfix = parse_param(options, 'savefile_postfix', []);

[Is,wl,bbl,A_gt_total,gt_names,endmembers,endmembers_r] = prepare_data(dataset);
% show_endmembers(endmembers, wl, gt_names, struct('bbl',bbl));
% show_ROIs_w_GT(Is,wl,select_image_indices,A_gt_total);

endmembers_ori = endmembers;
%% select good bands for spectra library
for i = 1:length(endmembers)
    endmembers{i} = endmembers{i}(:,bbl);
    endmembers_r{i} = endmembers_r{i}(:,bbl);
end

% endmembers = reduce_library(endmembers);

%% select specific images for running
if ~isempty(select_image_indices)
    Is = Is(select_image_indices);
    A_gt_total = A_gt_total(select_image_indices,:);
end
As = cell(1,length(Is));

%% transpose endmembers for BCM and AAM
switch algo
    case 'BCM'
        for j = 1:length(endmembers)
            endmembers{j} = max(min(endmembers{j}', 1), 0);
        end
    case 'AAM'
        for j = 1:length(endmembers_r)
            endmembers_r{j} = endmembers_r{j}';
        end
    case 'GMM'
        options.reduced_dim = parse_param(options,'reduced_dim',10);  
            
        % This is the only parameter that should be changed according
        % to the dataset. It controls the number of components in the
        % GMM distribution. It does not have a big impact on the
        % unmixing accuracy but has a big impact on the computational
        % effiency. Empirically, for 6 endmember classes or less, it is set to
        % (total number of spectra in the library) / 100000. However, 
        % it is better not to exceed 0.2. Hence, for too large
        % library, it is recommended to randomly discard spectra such
        % that the total number of spectra is less than 20000.
        options.thresh_CVIC = parse_param(options,'thresh_CVIC',0.05);          

        [mapping,components] = estimate_components_stable(endmembers, options);
        options.use_predefined_projection_components = 1;
        options.project_mapping = mapping;
        options.predefined_components = components;
%         show_scatter_and_wv_refl(Is,wl,endmembers_ori,bbl,gt_names,options);
    otherwise
end

%% batch processing all images

t_start = tic;

for i = 1:length(Is)    
    disp(['---- Start processing image No. ',num2str(i),' ----']);
    I = Is{i};
    I(:,:,~bbl) = [];
    Y = reshape_hsi(I);
    [rows,cols,B] = size(I);
    switch algo
        case 'GMM'
            D = 0.001^2 * eye(B);
            options.beta1 = parse_param(options,'beta1',0);
            options.beta2 = parse_param(options,'beta2',0);
            options.show_fig = 0;
            options.names = gt_names;
            options.D = D;
            
            [A,R,w_jk,mu_jk,sigma_jk,extra] = gmm_hu_ex(I, endmembers, options);
            % E = gmm_hu_endmember(I,A,D,w_jk,mu_jk,sigma_jk);
        case 'NCM' % GMM-1
            D = 0.001^2 * eye(B);
            options.beta1 = 0;
            options.beta2 = 0;
            options.show_fig = 0;
            options.names = gt_names;
            options.D = D;
            options.max_num_comp = 1;

            [A,R,w_jk,mu_jk,sigma_jk,extra] = gmm_hu_ex(I, endmembers, options);
        case 'NCM_sampling'
            num_endm = length(endmembers);
            E = zeros(B, num_endm);
            cov1 = zeros(B, B, num_endm);
            for j = 1:num_endm
                E(:,j) = mean(endmembers{j})';
                cov1(:,:,j) = cov(endmembers{j});
            end
            [Parameters] = unmixP_NCM_Parameters;
%             Parameters.NumberIterations = 3000;
            [P] = unmixP_NCM(Y',E,cov1,Parameters);
            A = P';
        case 'BCM'
            [Parameters] = BCMParameters(endmembers);

            disp('BCM-Spectral-QP Unmixing...');

            [A] = BCM(I, Parameters, 1);%BCM-Spectral-QP
        case 'AAM'
            [index1, A, reconstruction1, error1] = AAM(Y', endmembers_r);
            A = A';
        % This is a self-implemented MESMA. It is not the same as the
        % MESMA in the Viper tool (the one used in the paper).
        case 'MESMA'
            [A,E] = MESMA(I, endmembers_r);
        otherwise
    end
    As{i} = A;
    disp(['MAE of image No. ',num2str(i),' is ', ...
        num2str(mdiff(A_gt_total(i,:),mean(A,1),0))]);
%     close all;
end

t_elapsed = toc(t_start);
disp(['Total running time for ',algo,' is ', num2str(t_elapsed), ' seconds']);

if ischar(dataset)
    savefile = ['result_',dataset,'_',algo];
else
    savefile = ['result_custom_data_',algo];
end
save([savefile,savefile_postfix,'.mat'],'As','options');

error = calc_error_measure(A_gt_total, As);
total_time = t_elapsed;

end

