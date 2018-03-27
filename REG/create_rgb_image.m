function [rgb1, F, Sel, H] = create_rgb_image(I, wl, SRF, noise_rgb)
ideal_blue_wl = 470;
ideal_green_wl = 540;
ideal_red_wl = 650;

options = struct('ideal_red_wl',ideal_red_wl,'ideal_green_wl', ...
    ideal_green_wl,'ideal_blue_wl',ideal_blue_wl);
rgb_ind = find_rgb_ind(wl,options);
F = zeros(length(wl),3);

size_srf = floor(length(SRF) / 2);
for k = 1:size(F,2)
    start_ind = rgb_ind(k)-size_srf;
    end_ind = rgb_ind(k)+size_srf;
    inds = (start_ind : end_ind)';
    mask = inds >= 1 & inds <= size(F,1);
    start_ind = max(start_ind, 1);
    end_ind = min(end_ind, size(F,1));
    F(start_ind:end_ind, k) = SRF(mask);
end
Y = reshape(I, [size(I,1)*size(I,2) size(I,3)]);
rgb = Y * F;
rgb = reshape(rgb, [size(I,1),size(I,2),3]);

rgb1 = add_noise(rgb*255, noise_rgb*255);
rgb1(rgb1<0) = 0;
rgb1(rgb1>255) = 255;
rgb1 = uint8(rgb1);

F = F*255;
ind_sel = find(any(F~=0,2));
Sel = zeros(length(wl), length(ind_sel));
for i = 1:size(Sel,2)
    Sel(ind_sel(i),i) = 1;
end
H = F(ind_sel,:);
