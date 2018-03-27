function [h] = getHist(A,bitLen)

[N,dim] = size(A);
hashcode = zeros(N,1);
for d = 1:dim
    hashcode = hashcode*(2^bitLen) + A(:,d);
end

hashcode = sort(hashcode);
index = find((hashcode(1:N-1)-hashcode(2:N)) ~= 0);
h = ([index; N] - [0; index])/N;
