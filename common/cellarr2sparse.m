function A = cellarr2sparse(J2,A2,m,n)
%CELLARR2SPARSE Convert cell array of Js (cols) and Vs (value) to sparse
%matrix.
% Input:
%   J2 - cell array of Js (col) in which the index of J is the I (row)
%   A2 - cell array of values
%   m - rows of the sparse matrix
%   n - columns of the sparse matrix
% Output:
%   A - sparse matrix
total_size = 0;
for i = 1:length(J2)
    total_size = total_size + length(J2{i});
end

I1 = zeros(total_size,1);
J1 = zeros(total_size,1);
A1 = zeros(total_size,1);

start = 1;
for i = 1:length(J2)
    len = length(J2{i});
    if len == 0
        continue;
    end
    I1(start:start+len-1) = i * ones(len,1);
    J1(start:start+len-1) = J2{i}(:);
    A1(start:start+len-1) = A2{i}(:);
    start = start + len;
end
A = sparse(I1,J1,A1,m,n);
