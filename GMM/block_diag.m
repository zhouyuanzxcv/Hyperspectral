function Sigma1 = block_diag(Sigma, options)
%BLOCK_DIAG Construct a block diagonal sparse matrix.
% Input:
%   Sigma - B by B by N matrix where each page is placed as diagonal
%           element.
% Output:
%   Sigma1 - BN by BN sparse matrix
[B,~,N] = size(Sigma);

if nargin > 1
    Is = options.Is;
    Js = options.Js;
else
    Is = (1:N*B);
    Is = repmat(reshape(Is, [B,1,N]), 1, B);
    Is = Is(:);

    Js = (1:N*B);
    Js = repmat(reshape(Js, [1,B,N]), B, 1);
    Js = Js(:);
end

Sigma1 = sparse(Is, Js, Sigma(:), N*B, N*B);

end

