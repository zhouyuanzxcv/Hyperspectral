function [endmembers,I,Y,R_gt,A_gt,names,wl,Y_noise] = prepare_supervised_unmixing(dataset, options)
%PREPARE_SUPERVISED_UNMIXING Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    options = [];
end

if nargin < 1
    dataset = '03';
end

force_positive_I = 0;

if isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end

Y_noise = 0;

switch dataset
    case '00'
        load('toy_image_end_var_3.mat');
    case '001'
        load('toy_image_end_var_3.mat');
    case '002'
        load('toy_image_end_var_3_B.mat');
    case '01'
        if 0
            options = [];
            options.chosen = [1 3];
            [I,E_gt,A_gt,names,wl,Y_noise,extra] = generate_toy_image_unified(...
                60,60,'image','dirichlet',60,options);
            pure_spectra = extra.pure_spectra;
            save('toy_image_end_var_pavia_A_spectra_snr60.mat','I','E_gt',...
                'A_gt','names','wl','pure_spectra','Y_noise');
        else
%             load('toy_image_end_var_pavia_A_spectra_13_snr20.mat');
            load('C:\Data\GMM\Y_noise_10\toy_image_end_var_pavia_A_spectra_13_no_1.mat');
        end
    case '02'
        load('../GMM/toy_image_end_var_snr60_1231.mat');
    case '03'
        if 0
            options = [];
            options.chosen = [6 2];
            [I,E_gt,A_gt,names,wl,Y_noise,extra] = generate_toy_image_unified(...
                60,60,'gmm','dirichlet',60,options);
            ws = extra.ws;
            mus = extra.mus;
            sigmas = extra.sigmas;
            save('toy_image_end_var_dirichlet_snr60.mat','I','E_gt',...
                'A_gt','names','wl','ws','mus','sigmas','Y_noise');   
        else
            load('toy_image_end_var_dirichlet_snr60.mat');
        end
    case '21'
        load('../../Data/PaviaUniversity_A.mat');
        A_gt = double(A_gt);
    case '62'
        load('../../data/muufl_gulfport_B.mat');
        A_gt = double(A_gt);
    case '7'
        addpath('../../../HyperspectralImageAnalysis/Project2');
        [I,wl,pure_spectra,names] = load_santabar_spectra_lib('91');
        A_gt = zeros(size(I,1), size(I,2), length(pure_spectra));
        R_gt = zeros(length(pure_spectra), size(I,3));
    case '71'
        addpath('../../../HyperspectralImageAnalysis/Project2');
        [I,wl,pure_spectra,names] = load_santabar_spectra_lib('259');
        A_gt = zeros(size(I,1), size(I,2), length(pure_spectra));
        R_gt = zeros(length(pure_spectra), size(I,3));
    otherwise
end

if force_positive_I
    I(I<0) = 0;
    
    % for BCM, the first band of the gulfport dataset contains too many 0s
    % s.t. one endmember has all the pure pixels give 0 at the first band
    if strcmp(dataset,'62')
        I(:,:,1) = I(:,:,2);
        I(:,:,3) = I(:,:,4);
    end
end

[Y,A_gt,rows,cols] = reshape_hsi(I,A_gt);
[N,M] = size(A_gt);
[N,B] = size(Y);

switch dataset
    case {'0','01','7','71'} % pure_spectra is a cell array of endmember spectra
        endmembers = pure_spectra;
    case {'03'} % slice E_gt to get the endmember spectra
        endmembers = cell(1,M);
        for j = 1:M
            endmembers{j} = reshape(E_gt(j,:,:),B,N)';
        end
    otherwise % A_gt uses 1 to mark pure pixels
        endmembers = cell(1,M);
        for j = 1:M
            X = Y(A_gt(:,j)==1, :);
            endmembers{j} = X;
        end
end

R_gt = zeros(M,B);
for j = 1:M
    R_gt(j,:) = mean(endmembers{j});
end


