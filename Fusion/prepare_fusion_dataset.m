function [I_hyper,I_rgb,I_gt,g,G,extra,wl,wl_ori,Sel,s,bbl,H1] = ...
    prepare_fusion_dataset(dataset, options)
%PREPARE_FUSION_DATASET Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    options = [];
end

I_hyper = [];
I_rgb = [];

switch dataset
    case {'pavia_nonrigid_algo','pavia_nonrigid_gt','pavia_robustness'}
        if strcmp(dataset, 'pavia_nonrigid_gt')
            all_data = FusionDataLoader.PaviaFusion_GT_Reg();
        elseif strcmp(dataset, 'pavia_nonrigid_algo')
            all_data = FusionDataLoader.PaviaFusion_LSQfreeform_Reg();
        elseif strcmp(dataset, 'pavia_robustness')
            all_data = FusionDataLoader.PaviaFusion_robustness_LSQfreeform_Reg();
        end
        
        data_index = parse_param(options,'data_index',7);
        selected_data = all_data(data_index);
        
        I_gt = selected_data.I_gt;
        I_hyper = selected_data.I_hyper;
        wl = selected_data.wl;
        I_rgb = selected_data.I_rgb;
        s = selected_data.s; % s = (sx, sy)
        extra = selected_data.extra; % extra = [extra_x, extra_y]
        g = selected_data.g;
        G = selected_data.G;
        Sel = selected_data.S;
        H1 = selected_data.H1;        
        bbl = true(1,length(wl));
        
        wl_ori = wl;
    case 'gulfport'
        load('gulfport_fusion.mat');
        I_rgb = double(I_rgb);
        wl_ori = wl;
    case 'salton_sea'
        load('salton_sea_roi.mat');
        I1 = imread('salton_sea_roi_3.2015.jpg');
        load('reg_result_salton_sea.mat');
        bbl = logical(bbl);
        
        wl_ori = wl;
        [I_hyper,wl,I_rgb,g,H1,extra,s,G,Sel] = interpret_reg_params(I,wl,bbl,...
            I1,U2,V2,rho,s2,T2,sigma2);

        I_gt = zeros(size(I_rgb,1), size(I_rgb,2), length(bbl));
    otherwise
end

H1 = H1/255;
I_rgb = I_rgb/255;

% the competing methods do not have the extra margin boundary
if extra(1) > 0 
    I_rgb = remove_rgb_margin(I_rgb, extra);
    I_gt = remove_rgb_margin(I_gt, extra);
end

end

