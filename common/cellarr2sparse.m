function A = cellarr2sparse(J2,A2,m,n)
%CELLARR2SPARSE Summary of this function goes here
%   Detailed explanation goes here
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
