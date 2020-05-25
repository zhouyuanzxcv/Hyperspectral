function [ X ] = hyperFcls( M, U )
%HYPERFCLS Performs fully constrained least squares on pixels of M.
% hyperFcls performs fully constrained least squares of each pixel in M 
% using the endmember signatures of U.  Fully constrained least squares 
% is least squares with the abundance sum-to-one constraint (ASC) and the 
% abundance nonnegative constraint (ANC).
%
% Usage
%   [ X ] = hyperFcls( M, U )
% Inputs
%   M - HSI data matrix (p x N)
%   U - Matrix of endmembers (p x q)
% Outputs
%   X - Abundance maps (q x N)
%
% References
%   "Fully Constrained Least-Squares Based Linear Unmixing." Daniel Heinz, 
% Chein-I Chang, and Mark L.G. Althouse. IEEE. 1999.
%
% Written by
%    Matlab Hyperspectral Toolbox - Toolbox of advanced algorithms for
%    hyperspectral processing and exploitation.
%    http://sourceforge.net/apps/mediawiki/matlabhyperspec/index.php?title=Main_Page.
%    Toolbox Available at http://sourceforge.net/projects/matlabhyperspec/.

if (ndims(U) ~= 2)
    error('M must be a p x q matrix.');
end

[p1, N] = size(M);
[p2, q] = size(U);
if (p1 ~= p2)
    error('M and U must have the same number of spectral bands.');
end

p = p1;
X = zeros(q, N);
Mbckp = U;
for n1 = 1:N
    count = q;
    done = 0;
    ref = 1:q;
    r = M(:, n1);
    U = Mbckp;
    while not(done)
        als_hat = inv(U.'*U)*U.'*r;
        s = inv(U.'*U)*ones(count, 1);

        % IEEE Magazine method (http://www.planetary.brown.edu/pdfs/3096.pdf)
        % Contains correction to sign.  Error in original paper.
        afcls_hat = als_hat - inv(U.'*U)*ones(count, 1)*inv(ones(1, count)*inv(U.'*U)*ones(count, 1))*(ones(1, count)*als_hat-1);

        % See if all components are positive.  If so, then stop.
        if (sum(afcls_hat>0) == count)
            alpha = zeros(q, 1);
            alpha(ref) = afcls_hat;
            break;
        end
        % Multiply negative elements by their counterpart in the s vector.
        % Find largest abs(a_ij, s_ij) and remove entry from alpha.
        idx = find(afcls_hat<0);
        afcls_hat(idx) = afcls_hat(idx) ./ s(idx);
        [val, maxIdx] = max(abs(afcls_hat(idx)));
        maxIdx = idx(maxIdx);
        alpha(maxIdx) = 0;
        keep = setdiff(1:size(U, 2), maxIdx);
        U = U(:, keep);
        count = count - 1;
        ref = ref(keep);
    end
    X(:, n1) = alpha;
end

return;
