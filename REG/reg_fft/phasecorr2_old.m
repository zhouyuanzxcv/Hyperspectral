function d = phasecorr2(I1, I)
size_I1 = size(I1);
size_I = size(I);
outSize = size_I1(1:2) + size_I(1:2) - 1;
I = cat(3, ones(size(I,1), size(I,2)), I);
ds = [];

% [I,M] = scm_decompose(I);
% mask = create_gaussian(outSize(1), outSize(2));

for k = 1:size(I1,3)
    f = [];

    F1 = fft2(I1(:,:,k), outSize(1), outSize(2));
%     F1 = F1 .* mask;
        
    for i = 1:size(I,3)
        Fi = fft2(I(:,:,i), outSize(1), outSize(2));
%         Fi = Fi .* mask;
        
        Fi1 = Fi .* conj(F1) ./ (abs(F1+eps).^2);
        fi = ifft2(Fi1, 'symmetric');
        f = cat(3, f, fi);
    end
    
    A = reshape(f, [size(f,1)*size(f,2), size(f,3)]);
    [V,D] = eigs(A'*A, 1, 'lm');
    
    d = zeros(size(f,1), size(f,2));
    for i = 1:length(V)
        d = d + V(i) * f(:,:,i);
    end

    if max(-d(:)) > max(d(:))
        d = -d;
    end
    
    ds = cat(3, ds, d);
end

d = mean(ds, 3);

function [A,E] = scm_decompose(I)
I = I(:,:,2:end);
I = I/255;
options = [];
options.beta1 = 0.001;
options.beta2 = 0;
options.rho1 = 0.001;

options.show_figure = 0;
options.skip_uncertainty = 1;
options.use_single_variance_for_unmixing = 1;

M = 7;
[A,E,mu] = scm(I,M,options);
A = reshape(A, [size(I,1),size(I,2),M]);



function g = create_gaussian(rows, cols)
sigma_y = rows/4;
sigma_x = cols/4;
x = create_coordinate(cols);
y = create_coordinate(rows);
[X,Y] = meshgrid(x,y);
g = exp(-0.5*(X.^2/(sigma_x^2) + Y.^2/(sigma_y^2)));
g = g/sum(g(:));

function coord = create_coordinate(rows)
if mod(rows,2) == 1
    r = floor(rows/2);
    coord = (-r:r);
else
    r = rows/2;
    coord = (-r+0.5:r-0.5);
end