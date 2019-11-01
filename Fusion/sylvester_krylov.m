function X = sylvester_krylov(A, B, C, mode, options)
%SYLVESTER_KRYLOV Solve the sylvester equation AX + XB = C using the
%standard Krylov subspace method
% Input:
%   A - large sparse N by N matrix
%   B - small M by M matrix
%   C - N by M matrix
% Output:
%   X - N by M matrix
if nargin < 4
    mode = 'adaptive';
end
if nargin < 5
    options = [];
end

[U,Lambda] = eig(B);
lambdas = diag(Lambda);
D = C*U;
[N,N1] = size(A);
[M,M1] = size(B);
I = speye(N);
Y = zeros(N,M);
if strcmp(mode, 'original')
    for i = 1:M
        Y(:,i) = (A + lambdas(i)*I) \ D(:,i);
    end
elseif strcmp(mode, 'iterative') % when A is sparse but has too many bands, decomposing it takes memory
    for i = 1:M
        Y(:,i) = pcg(A + lambdas(i)*I, D(:,i));
    end
elseif strcmp(mode, 'adaptive') % when B has very few nonzero eigenvalues, this is much faster
    zero_inds = abs(lambdas) < 1e-15*max(abs(A(:)));
    % for many ds, decompose A by cholesky decomposition
    if isfield(options, 'factorized')
        Y(:,zero_inds) = options.factorized \ D(:,zero_inds);
    else
        Y(:,zero_inds) = A \ D(:,zero_inds);
    end
    for i = find(~zero_inds)'
%         Y(:,i) = (A + lambdas(i)*I) \ D(:,i);
        Y(:,i) = pcg(A + lambdas(i)*I, D(:,i));
    end
elseif strcmp(mode, 'adaptive_accurate')
    zero_inds = abs(lambdas) < 1e-15*max(abs(A(:)));
    % for many ds, decompose A by cholesky decomposition
    nonzero_inds = find(~zero_inds)';
    
    Y(:,zero_inds) = A \ D(:,zero_inds);
    for i = nonzero_inds
        Y(:,i) = (A + lambdas(i)*I) \ D(:,i);
    end
end
X = Y*U';
    
end

