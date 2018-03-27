function [ output_args ] = test_reg( input_args )
%TEST_REG Summary of this function goes here
%   Detailed explanation goes here
addpath('analysis');
addpath('reg_fft');
addpath('transform');
addpath('optimization');
addpath('mi_bspline');

close all

rho = 0;
lambda = 1e-3;

dataset = '20';
switch dataset
    case '00'
        noise_hyper = 0;
        noise_rgb = 0;
        [I1,rgb1,wl,T,s,sigma] = create_pavia_simulated();
        
        % save hyperspectral image and RGB image
        save('reg_pavia_dataset.mat','I1','wl','rgb1','T','s','sigma');
    case '10'
        load('reg_pavia_dataset.mat');
        U = zeros(size(I1,1),size(I1,2));
        V = U;
        options.s = 4.5;
        options.rho = 3;
    case '20'
        load('salton_sea_roi.mat');
%         rgb1 = imread('salton_sea_color.png');
        rgb1 = imread('salton_sea_roi_3.2015.jpg');
        bbl = logical(bbl);
        I1 = I(:,:,bbl);
        wl = wl(bbl);
        options.s = 10.4;
        % 20m is the IFOV of AVIRIS, 16.9m is the spatial resolution
        rho = ceil((20/16.9) * options.s / 2);
    case '30'
        load('neon_sjer_roi.mat');
        rgb1 = imread('2013_SJER_AOP_Camera_sample.tif');
        bbl = logical(bbl);
        I1 = I(:,:,bbl);
        wl = wl(bbl);
        options.s = 4;
    otherwise
end

options.rho = rho;
options.srf_range = 'color';
options.lambda = lambda;
options.show_figure = 1;

figure, imshow(uint8(rgb1));
figure, imshow(uint8(retrieve_rgb(I1, wl) * 255));

%% optimization
algo = 3;
switch algo
    case 0 % proposed nonrigid
        [T2,degree2,t2,s2,sigma2,U2,V2] = reg_hyper_rgb(rgb1, I1, wl, options);
    case 1 % proposed rigid
        options.reg_method = 'rigid';
        [T2,degree2,t2,s2,sigma2,U2,V2] = reg_hyper_rgb(rgb1, I1, wl, options);
    case 2 % MI rigid
        options.reg_metric = 'MI histogram';
        options.reg_method = 'rigid';
        [T2,degree2,t2,s2,sigma2,U2,V2] = reg_hyper_rgb(rgb1, I1, wl, options);
    case 3 % MI bspline
        options.reg_metric = 'MI histogram';
        options.nonrigid_model = 'MI bspline';
        [T2,degree2,t2,s2,sigma2,U2,V2] = reg_hyper_rgb(rgb1, I1, wl, options);
    otherwise
end
save('reg_result.mat','T2','degree2','t2','s2','sigma2','U2','V2','rho');

% resulting color image (scaled)
I2 = transform(double(rgb1), T2, s2, sigma2, size(I1,2), size(I1,1), options);
figure('name','final transformed color image');
imshow(uint8(I2));

% resulting hyperspectral image
options.isRigid = 0;
options.useTranslationField = 1;
options.U = U2;
options.V = V2;
I_hyper_reg = transform(I1,eye(3),[1,1],1e-3,size(I1,2),size(I1,1),options);
figure('name','final transformed hyperspectral image');
imshow(uint8(retrieve_rgb(I_hyper_reg, wl) * 255));


% calc_reg_error()


