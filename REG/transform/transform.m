function [I1,I2,g,extra] = transform(I, T, s, sigma, cols, rows, options)
% TRANSFORM
% Output:
%   I1 - final transformed image
%   I2 - image before scaling
%   g - point spread function
%   extra - [extra_x, extra_y] extra pixels in the PSF transform

% Coordinates of image are assumed to range from 0.5 to cols - 0.5.
if nargin < 7
    options = [];
    options.RaiseExceptionForExceeding = 0;
    options.FillValues = 0;
    options.isRigid = 1;
end

if isfield(options, 'rho') && options.rho > 0
    sigma = [sigma, options.rho];
end

s1 = ceil(s(1)); % scale_x
s2 = ceil(s(2)); % scale_y
cols1 = s1 * cols;
rows1 = s2 * rows;

delta_x = (s(1)*cols - 1) / (cols1 - 1);
delta_y = (s(2)*rows - 1) / (rows1 - 1);

[g,extra] = generate_PSF(sigma, s, delta_x, delta_y);

% do the (non)rigid transformation
viewport_x = [-extra(1)*delta_x+0.5, s(1)*cols+extra(1)*delta_x-0.5];
viewport_y = [-extra(2)*delta_y+0.5, s(2)*rows+extra(2)*delta_y-0.5];

if ~isfield(options,'isRigid') || options.isRigid % use the rigid transformation
    I2 = rigid_transform(I, T, viewport_x, viewport_y, cols1+2*extra(1), ...
        rows1+2*extra(2), options);
else  % use the nonrigid transformation
    if isfield(options, 'useGRBF') && options.useGRBF == 1
        options.centers = (options.centers - 0.5) .* s;
        options.sigmas = options.sigmas .* s;
    end
    I2 = nonrigid_transform(I, T, viewport_x, viewport_y, cols1+2*extra(1), ...
        rows1+2*extra(2), options);
end

% do the scaling with PSF
I1 = PSF_transform(I2, g, cols, rows, extra);


function [g,extra] = generate_PSF(sigma, s, delta_x, delta_y)
if length(sigma) == 2
    width = sigma(2);
    sigma = sigma(1);
    if width > 3*sigma
        width = 3*sigma;
    end
elseif isscalar(sigma)
    width = 2*sigma;
end
    
num = [ceil(s(1)), ceil(s(2))];
[x,extra_x] = calc_coord_PSF(num(1),delta_x,width);
[y,extra_y] = calc_coord_PSF(num(2),delta_y,width);

if length(x) == 1 && length(y) == 1
    g = 1;
else
    [X,Y] = meshgrid(x,y);
    f = exp(-(1/(2*sigma^2)) * (X.^2 + Y.^2));
    f(sqrt(X.^2 + Y.^2) > width) = 0;
    g = f./sum(f(:));
end
extra = [extra_x,extra_y];

function [coord,extra] = calc_coord_PSF(num, delta, width)
width = round(width / delta);
if width < num/2
    width = floor(num/2);
end

extra = width - floor(num/2);
if mod(num,2) == 1 % the window size is an odd number
    coord = (-width:width) * delta;
else 
    coord = (-width+0.5:width-0.5) * delta;
end

