function I1 = rigid_transform(I, T, viewport_x, viewport_y, x_num, y_num, options)
if nargin < 7
    options = [];
end
RaiseExceptionForExceeding = parse_param(options,'RaiseExceptionForExceeding',0);
FillValues = parse_param(options,'FillValues',0);

xi = linspace(viewport_x(1), viewport_x(2), x_num);
yi = linspace(viewport_y(1), viewport_y(2), y_num);
[Xi,Yi] = meshgrid(xi,yi);
I1 = rigid_transform_internal(I, T, Xi, Yi);

if RaiseExceptionForExceeding
    if any(isnan(I1(:)))
        throw(MException('Registration:ExceedBound', ...
            'Rigid transform samples points outside'));
    end
else
    I1(isnan(I1)) = FillValues;
end


function I1 = rigid_transform_internal(I, T, Xi, Yi)
%TRANSFORMIMAGE Transform image I(x) based on 3 by 3 transform matrix T on x.
% Generate coordinates in the rotated image
% Transform the coordinates
newCoord = T * [Xi(:)';Yi(:)';ones(size(Xi(:)'))];
newX = newCoord(1,:)./newCoord(3,:);
newY = newCoord(2,:)./newCoord(3,:);
XI = reshape(newX, size(Xi));
YI = reshape(newY, size(Yi));
% Interpolate the new coordinates
I1 = [];
x = linspace(0.5, size(I,2) - 0.5, size(I,2));
y = linspace(0.5, size(I,1) - 0.5, size(I,1));
[X,Y] = meshgrid(x,y);
for i = 1:size(I,3)
    I1 = cat(3, I1, interp2(X, Y, I(:,:,i), XI, YI, 'linear', NaN));
end

