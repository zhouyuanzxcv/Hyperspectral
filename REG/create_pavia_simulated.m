function [I1,rgb1,wl,T,s,sigma] = create_pavia_simulated(noise_hyper, ...
    noise_rgb,degree,U,V)
%CREATE_PAVIA_SIMULATED Summary of this function goes here
%   Detailed explanation goes here
rows = 80;
cols = 50;

if nargin < 1
    noise_hyper = 1e-3;
    noise_rgb = 1e-3;
    degree = 10;
    
    U = zeros(rows,cols);
    V = zeros(rows,cols);
end

load('../../data/PaviaUniversity_corrected.mat');

s = [4.4 4.5]'; % [x_scale, y_scale]
% sigma = 4.5/2;
% degree = 10, t = [96.74, 27.21]'
T = create_T(-degree, [-100 -10], 'inv'); % T is the transform of the coordinate

% I = double(I);

% create reference image
options = [];
options.rho = 3;

I = double(I);

sigma = 10;
I1 = transform(I, T, s, sigma, cols, rows, options);

options.isRigid = 0;
options.useTranslationField = 1;
options.U = -U;
options.V = -V;

I1 = transform(I1, eye(3), [1,1], 1e-2, cols, rows, options);


I1 = add_noise(I1, noise_hyper);
I1(I1<0) = 0;
I1(I1>1) = 1;

% figure, imshow(uint8(retrieve_rgb(I1, wl)*255));

% create rgb image
wl_width = mean(wl(2:end) - wl(1:end-1));
srf_size = floor((0.12 / wl_width) / 2); % SRF range is 120nm
SRF = fspecial('gaussian', [2*srf_size + 1, 1], srf_size / 3);
% SRF = calculate_Gaussian_srf(0.070, wl); % use fwhm = 70nm for each band
[rgb1, sel_mat] = create_rgb_image(I, wl, SRF, noise_rgb);

if 0 % plot SRF
    figure, hold on;
    plot(wl, sel_mat(:,1), 'r', 'linewidth', 1);
    plot(wl, sel_mat(:,2), 'g', 'linewidth', 1);
    plot(wl, sel_mat(:,3), 'b', 'linewidth', 1);
    legend('Red','Green','Blue');
end

rgb1 = rgb1(1:500,:,:);

end

function SRF = calculate_Gaussian_srf(fwhm, wl)
wl_width = mean(wl(2:end) - wl(1:end-1));
sigma = (fwhm/wl_width)/2.355; % for a Gaussian, fwhm = 2.355*sigma
SRF = fspecial('gaussian', [2*round(3*sigma) + 1, 1], sigma);
end