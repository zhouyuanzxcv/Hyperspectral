function sigmas = noise_est_mlr(I)
%NOISE_EST_L Summary of this function goes here
%   Detailed explanation goes here
if ndims(I) == 3
    [rows,cols,B] = size(I);
    Y = reshape(I, [rows*cols, B]);
else
    Y = I;
    [N,B] = size(Y);
end

sigmas = zeros(B,1);
epsilon = 1e-18;
bin_num = 30;
K = 2;
[A,J2,A2] = calc_adjacency_graph_cached(Y,K);
for k = 1:B
    y = Y(:,k);
    yplus = Y(:,min(k+1,B));
    yminus = Y(:,max(k-1,1));
    lsd = zeros(size(y));
    for i = 1:length(y)
        X = [1 yplus(i) yminus(i);ones(K,1) yplus(J2{i}) yminus(J2{i})];
        y1 = [y(i);y(J2{i})];
        b = regress(y1,X);
        lsd(i) = sqrt(sum((y1 - X*b).^2) / (K+1));
    end
%     hist(lsd);
    binrange = linspace(min(lsd)-epsilon, max(lsd)+epsilon, bin_num);
    bincount = histc(lsd,binrange);
    [~,ind] = max(bincount);
    sigmas(k) = (binrange(ind) + binrange(ind+1)) * 0.5;
end
    
