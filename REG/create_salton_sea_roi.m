function [ output_args ] = create_salton_sea_roi( input_args )
%CREATE_SALTON_SEA_ROI Summary of this function goes here
%   Detailed explanation goes here
filedir = 'C:\Data\Salton_Sea_f140331t01p00r13_refl\';
filepath = fullfile(filedir, 'f140331t01p00r13_corr_v1.hdr');
rgbpath = fullfile(filedir, 'f140331t01p00r13_sc01_RGB.jpeg');

[I_ori,wl,params] = read_cat_hdr(filepath);

bbl = params.bbl;
I_ori(:,:,~bbl) = NaN;

I_ori = I_ori / 10000;

rgb_ori = imread(rgbpath);

row_inds = 311:361;
col_inds = 628:683;
I = I_ori(row_inds,col_inds,:);
rgb = rgb_ori(row_inds,col_inds,:);

figure,imshow(uint8(retrieve_rgb(I,wl)*255));
figure,imshow(rgb);

wl = wl/1000;
save('salton_sea_roi.mat','I','wl','rgb','bbl','params');

end

