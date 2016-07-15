function [Y1,U,err] = reduce_dim(Y)
%REDUCE_DIM Summary of this function goes here
%   Detailed explanation goes here
errors = [];
for ndims = 1:20
    [Y1,U,err] = pca_dr(Y,ndims);
    errors = [errors;err];
    if ndims > 3 && errors(end-1)-errors(end) < 1e-3*(errors(1)-errors(2))
        break;
    end
end

function [Y1,U,err] = pca_dr(Y,ndims)
[mappedY, mapping] = pca(Y, ndims);
U = mapping.M;
Y1 = Y*U;
recon_Y = Y1*U';
err = mean(mean(abs(recon_Y - Y)));
