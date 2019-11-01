function I1 = remove_rgb_margin(I, extra)
I1 = I(extra(2)+1:end-extra(2), extra(1)+1:end-extra(1), :);
end
