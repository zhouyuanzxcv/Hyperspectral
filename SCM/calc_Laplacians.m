function [K,H] = calc_Laplacians(I,seg_map,M,eta,beta1,beta2)
[rows,cols,B] = size(I);
N = rows*cols;

if isempty(seg_map)
    [W,Neighbors] = image2graph(I,eta,1e-9);
else
    [W,Neighbors] = image2graph(double(seg_map),1e-3,1e-9);
end

D = diag(sum(W,2));
L = D - W;
L = sparse(L);
K = L - beta2/beta1*speye(N);

W = ones(M,M);
D = diag(sum(W,2));
H = D - W;

W = diag(ones(1,B-1),-1) + diag(ones(1,B-1),1);
D = diag(sum(W,2));
G = D - W;
G = sparse(G);

