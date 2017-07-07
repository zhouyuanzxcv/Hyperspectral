function [A,R] = scm_init(Y, M, init_mode, options)
%SCM_INIT Summary of this function goes here
%   Detailed explanation goes here
[N,B] = size(Y);

if init_mode == -1
    R = options.R_gt;
elseif init_mode == 0 % randomly pick pixels
    R = Y(round(rand(1,M)*(N-1) + 1),:);
elseif init_mode == 1 % pca
    R = zeros(M,B);
    [mappedX, mapping] = pca(Y,ceil(M/2));
    for i = 1:ceil(M/2)
        [~,ind] = sort(mappedX(:,i));
        R(i*2-1,:) = Y(ind(1),:);
        R(i*2,:) = Y(ind(end),:);
    end
    R = R(1:M,:);
elseif init_mode == 2 % k-means
    reduced_dim = 10;
    [Y1, mapping] = pca(Y, reduced_dim);
    options = [NaN,NaN,NaN,0]; % suppress output
    C = fcm(Y1,M,options);
    R = C*mapping.M' + repmat(mapping.mean,[M 1]);
elseif init_mode == 3
    R = vca(Y','endmembers',M);
    R = R';
end

A = Y*R'*inv(R*R'+eye(M)*1e-6);
A = project_to_simplex(A);
