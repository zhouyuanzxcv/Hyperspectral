function [QI, t_elapsed] = run_fusion_algo(dataset, algo, options)
%RUN_FUSION_ALGO Summary of this function goes here
%   Detailed explanation goes here
addpath('../REG');
addpath('../REG/transform');
addpath('../REG/optimization');
addpath('../REG/reg_fft');
addpath('../common');
addpath('../common/load_image');
addpath('analysis');

addpath(genpath('methods'));
addpath('Quality_Indices'); % from the pansharpening toolbox

if nargin < 3
    options = [];
end
if nargin < 1
%     dataset = 'pavia_nonrigid_gt';
%     dataset = 'pavia_nonrigid_algo';
%     dataset = 'pavia_robustness';

%     dataset = 'salton_sea';
    dataset = 'gulfport';
    
    algo = 'proposed';
%     algo = 'GFPCA';
%     algo = 'CNMF';
%     algo = 'Bayesian_Naive';
%     algo = 'Bayesian_Sparse';
%     algo = 'HySure';
end

[I_hyper,I_rgb,I_gt,g,G,extra,wl,wl_ori,Sel,s,bbl,H1] = ...
    prepare_fusion_dataset(dataset, options);

% for Bayesian fusion
ratio = s(1); % assume s(1) == s(2). Note that s = [s_x, s_y]
KerBlu = reshape(g, fliplr(s+2*extra));
start_pos = ceil([s(2)/2, s(1)/2]); % (start_r, start_c)
% overlap = find(sum(Sel, 2))';
F = Sel * H1(2:end,:);

curr_folder = fileparts(mfilename('fullpath'));


t_start = tic;
switch algo
    case 'proposed'
        
        % s = [s_x, s_y] is the scale. extra = [extra_x, extra_y] is the
        % extra pixel of the PSF excluding the scale, e.g. if the PSF
        % covers 9 x 9 pixels while s = [5, 5], then extra = [2, 2].
        options.s = s;
        options.extra = extra;

        % The following parameters correspond to the PSF and SRF.
        % If PSF and SRF are not available, simply setting has_PSF_SRF = 0
        % will work.         
        options.has_PSF_SRF = 1;
        options.g = g;
        options.H1 = H1;
        options.G = G;
        options.S = Sel;
        
        % The following parameters can be manually changed by uncommenting
        % them
%         options.gamma = 0.5; % balance between two data fidelity terms
%         options.beta = 1; % regularization parameter
%         options.K_neighbor = 3; % number of neighbors (K)
%         options.r2 = 15; % rho_2
%         options.r1 = 1; % rho_1
        
        I_recon = im_fusion(I_hyper,wl,I_rgb,options);
    case 'GFPCA'
        cd methods/GFPCA
        % manually tuned parameters for pavia LSQfreeform reg
        I_recon = GFPCA(I_hyper,I_rgb,9, size(KerBlu,1), 0.001^2);
    case 'CNMF'
        cd methods/CNMF
        I_recon = CNMF_fusion(I_hyper,I_rgb);
    case 'Bayesian_Naive'
        cd methods/BayesFusion
        setup;
        [I_recon]= BayesianFusion(I_hyper,I_rgb,F',KerBlu,ratio,'Gaussian',start_pos);
    case 'Bayesian_Sparse'
        cd methods/BayesFusion
        setup;
        [I_recon]= BayesianFusion(I_hyper,I_rgb,F',KerBlu,ratio,'Sparse',start_pos);
    case 'HySure'
        cd methods/HySure
        I_recon = HySure_wrapper(I_hyper, I_rgb, KerBlu, ratio, F', start_pos);
    otherwise
end
t_elapsed = toc(t_start);
disp(['Time cost for ',algo,' is ',num2str(t_elapsed)]);

cd(curr_folder);

% Add bad bands
I_recon1 = nan(size(I_recon,1), size(I_recon,2), length(bbl));
I_recon1(:,:,bbl) = I_recon;
I_recon = I_recon1;

QI = QualityIndices(I_recon,I_gt,ratio);

if strcmp(dataset, 'salton_sea') || strcmp(dataset,'gulfport') % show fusion result
    savefile = ['result_',dataset,'_',algo,'.mat'];
    save_fusion_for_viewing(I_recon,wl_ori,bbl,savefile);
    analyze_hsi('UserData',savefile);
end

end

% function [Sel,H1,g,G] = change_interface(s, extra, KerBlu, F, I_hyper, I_rgb)
% %If PSF and SRF are available but provided as KerBlu
% %and F, the following code can convert them to H1 and G. In this
% %case, the start_pos should be start_pos = ceil([s(2)/2, s(1)/2]).
% [rows,cols,B] = size(I_hyper);
% [rows1,cols1,b] = size(I_rgb);
% [r_g,c_g] = size(KerBlu);
% assert(all(extra == [c_g - s(1),r_g - s(2)]/2));
% Sel = eye(B);
% H1 = [zeros(1,b);F];
% g = KerBlu(:);
% Cs = construct_C(I_hyper,I_rgb,s,extra);
% C = cat(2,Cs{:});
% G = g2G(C,g,rows*cols);
% 
% end


