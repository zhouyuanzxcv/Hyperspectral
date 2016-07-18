function [W,M] = image2graph(im,eta,threshold,fcn_weight)
%LABEL2GRAPH Construct a graph from the image based on one order
%neighborhood.
% Input: 
%   im - rows by cols by B image
%   eta - used in the heat kernel for calculating the weights
%   threshold - truncate small values in the weights to 0
%   fcn_weight - function handle to calculate the weight, replace the heat
%                kernel if inputted.
% Output:
%   W - weighted neighborhood matrix
%   M - indicator matrix for neighborhood
if nargin < 3
    threshold = 1e-9;
end

[m,n,B] = size(im);
eta = eta*sqrt(B);
N = m*n;

if nargin < 4
    fcn_weight = @(y1,y2) calc_weight(y1,y2,eta);
end

[W,M] = construct_W1(im,fcn_weight);

[I,J,A] = find(W);

T = max(A(:)) * threshold;
A1 = A>T;
I = I(A1);
J = J(A1);
A = A(A1);

W = sparse(I,J,A,N,N);

function [W,M] = construct_W1(im,fcn_weight)
[m,n,B] = size(im);
N = m*n;

[X,Y] = meshgrid(1:n,1:m);
X = X(:);
Y = Y(:);

I = [];
J = [];

X_top = X;
Y_top = Y - 1;
ind_invalid = Y_top < 1;
[I1,J1] = sub2ind1(X,Y,X_top,Y_top,ind_invalid,m,n);
I = [I;I1];
J = [J;J1];

X_bottom = X;
Y_bottom = Y + 1;
ind_invalid = Y_bottom > m;
[I1,J1] = sub2ind1(X,Y,X_bottom,Y_bottom,ind_invalid,m,n);
I = [I;I1];
J = [J;J1];

X_left = X - 1;
Y_left = Y;
ind_invalid = X_left < 1;
[I1,J1] = sub2ind1(X,Y,X_left,Y_left,ind_invalid,m,n);
I = [I;I1];
J = [J;J1];

X_right = X + 1;
Y_right = Y;
ind_invalid = X_right > n;
[I1,J1] = sub2ind1(X,Y,X_right,Y_right,ind_invalid,m,n);
I = [I;I1];
J = [J;J1];

im1 = reshape(im, m*n, B);

% A = calc_weight(im1(I,:)',im1(J,:)',eta);
A = fcn_weight(im1(I,:),im1(J,:));

W = sparse(I,J,A,N,N);
M = sparse(I,J,ones(size(I)),N,N);

function [I1,J1] = sub2ind1(X,Y,X1,Y1,ind_invalid,m,n)
X1(ind_invalid) = [];
Y1(ind_invalid) = [];
X(ind_invalid) = [];
Y(ind_invalid) = [];

I1 = sub2ind([m n],Y,X);
J1 = sub2ind([m n],Y1,X1);

function w = calc_weight(y1,y2,eta)
w = exp(-sum((y1 - y2).^2, 2) / (2*eta^2));

% obsolete
function W = construct_W(im,eta)
[m,n,B] = size(im);
N = m*n;

max_size = N*4;
I = zeros(max_size,1);
J = zeros(max_size,1);
A = zeros(max_size,1);
I_ind = 1;
J_ind = 1;
A_ind = 1;

for i = 1:m
    for j = 1:n
        idx1 = sub2ind([m,n],i,j);
        L1 = im(i,j,:);
        if i-1 >= 1
            I(I_ind) = idx1;
            J(J_ind) = sub2ind([m,n],i-1,j);
            A(A_ind) = calc_weight(L1,im(i-1,j,:),eta);
            I_ind = I_ind + 1;
            J_ind = J_ind + 1;
            A_ind = A_ind + 1;
        end
        if i+1 <= m
            I(I_ind) = idx1;
            J(J_ind) = sub2ind([m,n],i+1,j);
            A(A_ind) = calc_weight(L1,im(i+1,j,:),eta);
            I_ind = I_ind + 1;
            J_ind = J_ind + 1;
            A_ind = A_ind + 1;
        end
        if j-1 >= 1
            I(I_ind) = idx1;
            J(J_ind) = sub2ind([m,n],i,j-1);
            A(A_ind) = calc_weight(L1,im(i,j-1,:),eta);
            I_ind = I_ind + 1;
            J_ind = J_ind + 1;
            A_ind = A_ind + 1;
        end
        if j+1 <= n
            I(I_ind) = idx1;
            J(J_ind) = sub2ind([m,n],i,j+1);
            A(A_ind) = calc_weight(L1,im(i,j+1,:),eta);
            I_ind = I_ind + 1;
            J_ind = J_ind + 1;
            A_ind = A_ind + 1;
        end
    end
end

I = I(1:I_ind-1);
J = J(1:J_ind-1);
A = A(1:A_ind-1);
W = sparse(I,J,A,N,N);
