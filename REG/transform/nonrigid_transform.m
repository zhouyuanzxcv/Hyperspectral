function I1 = nonrigid_transform(I, T, viewport_x, viewport_y, x_num, y_num, options)
if nargin < 7
    options = [];
    options.centers = [];
    options.weights = [];
    options.useGRBF = 0;
    options.useTranslationField = 1;
    options.U = zeros(y_num, x_num);
    options.V = zeros(y_num, x_num);
end
RaiseExceptionForExceeding = parse_param(options,'RaiseExceptionForExceeding',1);
FillValues = parse_param(options,'FillValues',0);

xi = linspace(viewport_x(1), viewport_x(2), x_num);
yi = linspace(viewport_y(1), viewport_y(2), y_num);
[Xi,Yi] = meshgrid(xi,yi);

if isfield(options,'useGRBF') && options.useGRBF == 1
    I1 = nonrigid_transform_grbf(I, T, Xi, Yi, options.centers, options.weights);
elseif isfield(options,'useTranslationField') && options.useTranslationField == 1
    I1 = nonrigid_transform_translation_field(I, T, Xi, Yi, options.U, options.V);
end

if RaiseExceptionForExceeding
    if any(isnan(I1(:)))
        throw(MException('Registration:ExceedBound', ...
            'Rigid transform samples points outside'));
    end
else
    I1(isnan(I1)) = FillValues;
end


function I1 = nonrigid_transform_translation_field(I, T, Xi, Yi, U, V)
[XI, YI] = rigid_coord_transform(T, Xi, Yi);
XI = XI + U;
YI = YI + V;

xmin = 0.5;
xmax = size(I,2) - 0.5;
ymin = 0.5;
ymax = size(I,1) - 0.5;

x = linspace(xmin, xmax, size(I,2));
y = linspace(ymin, ymax, size(I,1));
[X,Y] = meshgrid(x,y);

I1 = interp2_hyper(I,X,Y,XI,YI,[xmin,xmax],[ymin,ymax]);


function I1 = nonrigid_transform_grbf(I, T, Xi, Yi, centers, weights)
%TRANSFORMIMAGE Transform image I(x) based on 3 by 3 transform matrix T on x.
range = round(3*sigma);
x = (-range:range);
y = (-range:range);
[X,Y] = meshgrid(x,y); 
grbf_mask = exp(-(X.^2 + Y.^2)/(2*sigma^2));

% Transform the coordinates rigidly
[XI, YI] = rigid_coord_transform(T, Xi, Yi);

% Transform the coordinates by GRBF
[cx, cy] = rigid_coord_transform(T, centers(:,1), centers(:,2));
for i = 1:size(centers,1)
    x1 = cx(i) + x;
    y1 = cy(i) + y;
    
    logical_ind = x1 > 0 & x1 <= size(XI,2);
    x1 = x1(logical_ind);
    grbf_mask1 = grbf_mask(:,logical_ind);
    
    logical_ind = y1 > 0 & y1 <= size(XI,1);
    y1 = y1(logical_ind);
    grbf_mask1 = grbf_mask1(logical_ind,:);
    
    XI(y1,x1) = XI(y1,x1) + grbf_mask1 * weights(i,1);
    YI(y1,x1) = YI(y1,x1) + grbf_mask1 * weights(i,2);
end
    
% Interpolate the new coordinates
I1 = [];
x = linspace(0.5, size(I,2) - 0.5, size(I,2));
y = linspace(0.5, size(I,1) - 0.5, size(I,1));
[X,Y] = meshgrid(x,y);
for i = 1:size(I,3)
    I1 = cat(3, I1, interp2(X, Y, I(:,:,i), XI, YI, 'linear', NaN));
end


function [XI, YI] = rigid_coord_transform(T, Xi, Yi)
newCoord = T * [Xi(:)';Yi(:)';ones(size(Xi(:)'))];
newX = newCoord(1,:)./newCoord(3,:);
newY = newCoord(2,:)./newCoord(3,:);
XI = reshape(newX, size(Xi));
YI = reshape(newY, size(Yi));
