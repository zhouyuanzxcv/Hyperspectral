function negative_MI = calc_mi_kdensity(img1,img2,index)
%CALC_MI_KDENSITY Summary of this function goes here
%   Detailed explanation goes here
n = sqrt(size(img1,1)*size(img1,2));
data = [img1(:),img2(:)];
[bandwidth,density,X,Y]=kde2d(data,n);
h1 = sum(density,1);
h2 = sum(density,2);
h12 = density(:);

e1 = -sum(h1.*log(h1));
e2 = -sum(h2.*log(h2));

e12 = -sum(h12.*log(h12));

negative_MI = e12 - e1- e2; % negative MI

end

