function C = construct_C(I,rgb,s,extra)
% Contruct a cell array of C_n correspondence matrices.
% Input:
%   I - rows by cols by B hyperspectral image
%   rgb - rows1 by cols1 by b color image
%   s - assumed to be a vector of 2 integers, (sx, sy)
%   extra - (extra_x, extra_y)
assert(extra(1)*2 + size(I,2)*s(1) == size(rgb,2));
assert(extra(2)*2 + size(I,1)*s(2) == size(rgb,1));

s1 = s + 2*extra;
R = prod(s1);

[rows,cols,B] = size(I);
[rows1,cols1,b] = size(rgb);
N = rows*cols;

C = cell(1,N);
N1 = rows1*cols1;

for i = 1:N
    I = (1:R);
    [r,c] = ind2sub([rows,cols],i);
    x1_inds = (c-1)*s(1)+1:c*s(1)+2*extra(1);
    y1_inds = (r-1)*s(2)+1:r*s(2)+2*extra(2);
    [X,Y] = meshgrid(x1_inds, y1_inds);
    J = sub2ind([rows1,cols1],Y(:),X(:))';
    V = ones(1,R);
    C{i} = sparse(I,J,V,R,N1);
end
