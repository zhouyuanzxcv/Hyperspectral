function [ output_args ] = test_fft_reg( input_args )
%TEST_FFT_REG Summary of this function goes here
%   Detailed explanation goes here
addpath('reg_fft');

close all

%% load or create image
dataset = '11';
switch dataset
    case '00'
        load('../../data/PaviaUniversity_corrected.mat');
        SRF = [0.2,0.4,0.6,0.4,0.2];
        [rgb,sel_mat] = create_rgb_image(I, wl, SRF);
        I1 = double(rgb);
        
        I = I(30:250,80:250,:);
        save('reg_pavia_rot_tran.mat','s','lambda','T','I','wl','I1');
    case '01'
        load('../../data/PaviaUniversity_corrected.mat');
        SRF = [0.2,0.4,0.6,0.4,0.2];
        [rgb,sel_mat] = create_rgb_image(I, wl, SRF);
        I1 = double(rgb);
        
        rows = 500;
        cols = 240;
        s = [1 1]'; % [x_scale, y_scale]
        lambda = 0.01;
        T = create_T(-10, [-100 0], 'inv'); % T is the transform of the coordinate

        % create reference image
        I = double(I);
        I = transform(I, T, s, lambda, cols, rows);
        I = I;
        save('reg_pavia_rot_tran.mat','s','lambda','T','I','wl','I1');
    case '02'
        load('../../data/PaviaUniversity_corrected.mat');
        SRF = [0.2,0.4,0.6,0.4,0.2];
        [rgb,sel_mat] = create_rgb_image(I, wl, SRF);
        I1 = double(rgb(:,:,1));
        
        rows = 500;
        cols = 240;
        s = [1 1]'; % [x_scale, y_scale]
        lambda = 0.01;
        T = create_T(-10, [-100 0], 'inv'); % T is the transform of the coordinate

        % create reference image
        I = transform(I1, T, s, lambda, cols, rows);
        I = I / 255;
        save('reg_pavia_rot_tran_1band.mat','s','lambda','T','I','wl','I1');
    case '11'
        load('reg_pavia_rot_tran.mat');
    case '12'
        load('reg_pavia_rot_tran_1band.mat');
    otherwise
end


I = I * 255;

% show test image
figure,imshow(uint8(I1));
% show reference image
if size(I,3) > 3
    figure,imshow(uint8(retrieve_rgb(I,wl)*255));
else
    figure,imshow(uint8(I));
end

%% Do phase correlation
% reg_imregcorr(I1, I, wl);
% d = findTranslationPhaseCorr(I1, I);
[angle, t, s] = reg_hyper_fft(I1, I, wl, 'rigid');
T = create_T(angle, t);
options = [];
options.RaiseExceptionForExceeding = 0;
I2 = transform(I1, T, s, lambda, size(I,2), size(I,1), options);
figure,imshow(uint8(I2));

function reg_imregcorr(I1, I, wl)
if size(I,3) > 3
    rgb_ind = find_rgb_ind(wl);
else
    rgb_ind = 1;
end

tform = imregcorr(I1(:,:,1), I(:,:,rgb_ind(1)), 'rigid','window',1);
I2 = imwarp(I1, tform, 'OutputView', imref2d(size(I)));
figure,imshow(uint8(I2));

function ref = create_ref_image(I)
T = create_T(0, [100;100], 'inv');
[rows,cols,~] = size(I);
try
    ref = rigid_transform(I, T, cols, rows, cols, rows);
catch err
    if strcmp(err.identifier, 'Registration:ExceedBound')
        disp('Set outside values to zeros');
    else
        throw(err);
    end
end



