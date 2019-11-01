function I_rgb1 = pad_rgb_margin(I_rgb, extra)
%PAD_RGB_MARGIN Replicate boundary values to pad image by extra size
I_rgb1 = zeros(size(I_rgb,1)+2*extra(2), size(I_rgb,2)+2*extra(1), size(I_rgb,3));
I_rgb1(extra(2)+1:end-extra(2), extra(1)+1:end-extra(1), :) = I_rgb;
I_rgb1(extra(2)+1:end-extra(2), 1:extra(1), :) = repmat(I_rgb(:,1,:),[1 extra(1) 1]);
I_rgb1(extra(2)+1:end-extra(2), end-extra(1)+1:end, :) = repmat(I_rgb(:,end,:), [1 extra(1) 1]);
I_rgb1(1:extra(2), :, :) = repmat(I_rgb1(extra(2)+1,:,:), [extra(2) 1 1]);
I_rgb1(end-extra(2)+1:end,:,:) = repmat(I_rgb1(end-extra(2),:,:), [extra(2) 1 1]);


end

