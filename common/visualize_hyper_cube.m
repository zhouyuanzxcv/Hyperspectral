function visualize_hyper_cube(I, wl, scale)
%VISUAL_HYPER_IMAGE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 1
    if 0
        I = imread('salton_sea_roi_3.2015.jpg');
        I = im2double(I);
        I = cat(3, I(:,:,3), I(:,:,2), I(:,:,1));
        wl = [0.47, 0.56, 0.65];
        scale = 1;
    else
        load('salton_sea_roi.mat');
        scale = 100 / size(I,1);
%         I = I;
    end
end

I1 = [];
for k = 1:size(I,3)
    I1(:,:,k) = imresize(I(:,:,k), scale, 'nearest');
end

rgb1 = retrieve_rgb(I1,wl);
[rows,cols,B] = size(I1);

rows1 = rows + round(rows*0.5);
cols1 = cols + round(rows*0.5);
rows_end = rows1 - rows + 1;
cols_end = cols1 - cols + 1;

rgb_full = ones(rows1,cols1,3);
select_pos_idx = linspace(1,rows_end,B+1);
for i = 1:length(select_pos_idx)-1
    if mod(i,20) == 0
        i
    end
    row_i = round(select_pos_idx(i));
    col_j = round(select_pos_idx(i));
    I2 = zeros(size(I1));
    inds = B-i+1;
%     inds(inds < 1 | inds > B) = [];
    I2(:,:,inds) = I1(:,:,inds);
    rgb2 = hyper2rgb(I2,wl);
    rgb_full(row_i:row_i+rows-1, col_j:col_j+cols-1, :) = rgb2;
end
rgb_full(rows_end:rows_end+rows-1, cols_end:cols_end+cols-1, :) = rgb1;

figure,imshow(rgb_full);

% imwrite(rgb_full, 'tmp.jpg');
end

function rgb = hyper2rgb(I1, wl)
[rows,cols,B] = size(I1);

options = [];
if wl(end) > 2 % more than 2 micrometer
    options.ideal_red_wl = 800;
    options.ideal_green_wl = 650;
    options.ideal_blue_wl = 470;
end

rgb_ind = find_rgb_ind(wl, options);
wl_delta = mean(wl(2:end) - wl(1:end-1)); 

if wl(end) > 2
    % use 400 num SRF
    half_step = round(0.2 / wl_delta);
else
    half_step = round(0.05 / wl_delta);
end
gauss = fspecial('gaussian', [2*half_step+1, 1], half_step/2);
gauss = gauss / max(gauss);
if wl(end) > 2
    gauss = 2 * gauss;
end

rgb_inds = min(rgb_ind) - half_step : max(rgb_ind) + half_step;
rgb_inds(rgb_inds < 1) = [];

% use SRF to visualize the color bands
rgb = zeros(size(I1,1),size(I1,2),3);
for k = 1:length(rgb_ind)
    inds = (rgb_ind(k) - half_step: rgb_ind(k) + half_step);
    inds_omit = inds < 1 | inds > size(I1,3);
    inds(inds_omit) = [];
    gauss1 = gauss;
    gauss1(inds_omit) = [];
    I1_1 = reshape(I1(:,:,inds), rows*cols, length(inds)) * gauss1;
    rgb(:,:,k) = reshape(I1_1, rows, cols);
end

rgb = rgb;

% if beyond visible bands
if mean(rgb(:)) < 0.1
    I2 = sum(I1,3);
    rgb = rgb + cat(3,I2,0.8*I2,0.8*I2);
end


end