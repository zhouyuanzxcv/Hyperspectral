function [I,R_gt,wl] = remove_noisy_bands(I,R_gt,wl,remove_bands)
%REMOVE_ Summary of this function goes here
%   Detailed explanation goes here
I(:,:,remove_bands) = [];
R_gt(:,remove_bands) = [];
wl(remove_bands) = [];

end

