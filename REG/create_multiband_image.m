function [I1, F, Sel, H] = create_multiband_image(I, wl, SRF_shape, SRF_centers, noise)
b = length(SRF_centers);
for i = 1:b
    [~,ind] = min(abs(wl - SRF_centers(i)));
    mb_inds(i) = ind;
end

F = zeros(length(wl),b);

if mod(length(SRF_shape),2) == 0 % if even number, make it odd
    SRF_shape = cat(1,SRF_shape,0);
end
    
size_srf = floor(length(SRF_shape) / 2);
for k = 1:size(F,2)
    start_ind = mb_inds(k)-size_srf;
    end_ind = mb_inds(k)+size_srf;
    inds = (start_ind : end_ind)';
    mask = inds >= 1 & inds <= size(F,1);
    F(inds(mask), k) = SRF_shape(mask);
end
Y = reshape(I, [size(I,1)*size(I,2) size(I,3)]);
Y1 = Y * F;
I1 = reshape(Y1, [size(I,1),size(I,2),b]);

I1 = add_noise(I1, noise);


ind_sel = find(any(F~=0,2));
Sel = zeros(length(wl), length(ind_sel));
for i = 1:size(Sel,2)
    Sel(ind_sel(i),i) = 1;
end
H = F(ind_sel,:);

end

