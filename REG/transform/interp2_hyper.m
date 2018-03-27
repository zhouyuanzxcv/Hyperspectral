function I1 = interp2_hyper(I,X,Y,XI,YI,xrange,yrange)
%INTERP2_HYPER Summary of this function goes here
%   Detailed explanation goes here
% xrange = [min(X(:)), max(X(:))]
% yrange = [min(Y(:)), max(Y(:))]

I1 = zeros(size(XI,1),size(XI,2),size(I,3));

% use replicate for the boundary condition of I
XI(XI < xrange(1)) = xrange(1);
XI(XI > xrange(2)) = xrange(2);
YI(YI < yrange(1)) = yrange(1);
YI(YI > yrange(2)) = yrange(2);

% this part is used for nonrigid transformation of hyperspectral images.
% For the time being, the interp2 is set to spline.
for i = 1:size(I,3)
    I1(:,:,i) = interp2(X, Y, I(:,:,i), XI, YI, 'spline');
end


end

