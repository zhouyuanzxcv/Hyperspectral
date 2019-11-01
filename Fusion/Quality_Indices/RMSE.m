function out = RMSE(ref,tar)
%--------------------------------------------------------------------------
% Root mean squared error (RMSE)
%
% USAGE
%   out = RMSE(ref,tar)
%
% INPUT
%   ref : reference HS data (rows,cols,bands)
%   tar : target HS data (rows,cols,bands)
%
% OUTPUT
%   out : RMSE (scalar)
%
%--------------------------------------------------------------------------
[rows,cols,bands] = size(ref);
out = (sum(sum(sum((tar-ref).^2)))/(rows*cols*bands)).^0.5;