function RGB = retrieve_rgb(data, wl, ideal_red_wl, ideal_green_wl, ideal_blue_wl)
%RETRIEVE_RGB Summary of this function goes here
%   Detailed explanation goes here

if mean(wl) < 10 % the unit is micrometer
    wl = wl*1e3;
end

if nargin < 3
    ideal_blue_wl = 470;
    ideal_green_wl = 540;
    ideal_red_wl = 650;
end

blue_ind = find(wl > ideal_blue_wl, 1, 'first');
green_ind = find(wl > ideal_green_wl, 1, 'first');
red_ind = find(wl > ideal_red_wl, 1, 'first');

blue_wl = wl(blue_ind);
green_wl = wl(green_ind);
red_wl = wl(red_ind);

disp(['Use ',num2str(blue_wl),'nm for blue, ', ...
    num2str(green_wl),'nm for green, ',num2str(red_wl),'nm for red.']);

I = data;
R = I(:,:,red_ind);
G = I(:,:,green_ind);
B = I(:,:,blue_ind);

RGB = cat(3,R,G,B);

RGB = RGB / max(RGB(:));


