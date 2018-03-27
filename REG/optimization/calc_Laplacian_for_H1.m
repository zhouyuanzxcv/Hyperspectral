function L1 = calc_Laplacian_for_H1(B1)
B = B1 - 1;
L = calc_Laplacian_for_H(B);
L1 = zeros(B1,B1);
L1(2:end,2:end) = L;

function L = calc_Laplacian_for_H(B)
W = diag(ones(1,B-1),-1) + diag(ones(1,B-1),1);
D = diag(sum(W,2));
L = D - W;
