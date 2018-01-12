function RGB = retrieve_rgb(data, wl, ideal_red_wl, ideal_green_wl, ideal_blue_wl, verbose)
%RETRIEVE_RGB Summary of this function goes here
%   Detailed explanation goes here
if nargin < 6
    verbose = 0;
end

if mean(wl) < 10 % the unit is micrometer
    wl = wl*1e3;
end

if nargin < 3
    ideal_blue_wl = 470;
    ideal_green_wl = 540;
    ideal_red_wl = 650;
end

options = struct('ideal_red_wl',ideal_red_wl,'ideal_green_wl', ...
    ideal_green_wl,'ideal_blue_wl',ideal_blue_wl);
rgb_ind = find_rgb_ind(wl,options);

red_ind = rgb_ind(1);
green_ind = rgb_ind(2);
blue_ind = rgb_ind(3);

blue_wl = wl(blue_ind);
green_wl = wl(green_ind);
red_wl = wl(red_ind);

if verbose
    disp(['Use ',num2str(blue_wl),'nm for blue, ', ...
        num2str(green_wl),'nm for green, ',num2str(red_wl),'nm for red.']);
end

I = data;
R = I(:,:,red_ind);
G = I(:,:,green_ind);
B = I(:,:,blue_ind);

RGB = cat(3,R,G,B);

RGB = RGB / max(RGB(:));


