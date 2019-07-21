function [rgb1, F, Sel, H] = create_rgb_image(I, wl, SRF_shape, noise_rgb)
ideal_blue_wl = 470;
ideal_green_wl = 540;
ideal_red_wl = 650;

rgb_centers = [ideal_red_wl, ideal_green_wl, ideal_blue_wl] / 1000;

[rgb, F, Sel, H] = create_multiband_image(I, wl, SRF_shape, rgb_centers, noise_rgb);

rgb1 = rgb*255;
rgb1(rgb1<0) = 0;
rgb1(rgb1>255) = 255;
rgb1 = uint8(rgb1);

F = F*255;
H = H*255;
end

