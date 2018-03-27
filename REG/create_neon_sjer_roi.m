function [ output_args ] = create_salton_sea_roi( input_args )
%CREATE_SALTON_SEA_ROI Summary of this function goes here
%   Detailed explanation goes here
filedir = 'C:\Data\NEON_AOP_sample_data_v2';
filepath = fullfile(filedir, 'Spectrometer\2013_SJER_AOP_NIS_sample.hdr');
% rgbpath = fullfile(filedir, 'RGB_Camera\2013_SJER_AOP_Camera_sample.tif');

[I_ori,wl,params] = read_cat_hdr(filepath);

% bbl = params.bbl;
% I_ori(:,:,~bbl) = NaN;
I_ori(I_ori == params.data_ignore_value) = NaN;
bbl = squeeze(any(any(isnan(I_ori),1),2));
bbl = logical(~bbl);

I_ori = I_ori / 10000;

% rgb_ori = imread(rgbpath);

row_inds = 1:250;
col_inds = 1:250;
I = I_ori(row_inds,col_inds,:);
% rgb = rgb_ori(row_inds,col_inds,:);

figure,imshow(uint8(retrieve_rgb(I,wl)*255));
% figure,imshow(rgb);

% wl = wl/1000;
save('neon_sjer_roi.mat','I','wl','bbl','params');

end

