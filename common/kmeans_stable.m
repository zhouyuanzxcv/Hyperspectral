function C = kmeans_stable(Y, M, num_kmeans)
% Find the most frequents centers given by kmeans.
%   Input:
%       Y - N by D data
%       M - number of clusters
%       num_kmeans - number of runs
%   Output:
%       C - Centers of the most frequent kmeans
[N,B] = size(Y);

%% reduce the dimension
reduced_dim = 10;
use_reduced_dim = 0;
if B > reduced_dim
    use_reduced_dim = 1;
    [Y, mapping] = pca(Y, reduced_dim);
    B = reduced_dim;
end


%% run multiple times
Cs = zeros(M,B,num_kmeans);
err_inds = [];
for i = 1:num_kmeans
    try
        Y1 = Y;
        options = [NaN,NaN,NaN,0]; % suppress output
%         C = fcm(Y1,M,options);
        [~, C] = kmeans(Y1, M, 'start', 'sample');
        Cs(:,:,i) = C;
    catch me
        err_inds = [err_inds,i];
    end
end
Cs(:,:,err_inds) = [];

C = find_most_frequent(Cs);

if use_reduced_dim
    C = C*mapping.M' + repmat(mapping.mean,[M 1]);
end

function C = find_most_frequent(Cs)
reduced_dim = 3;

[M,B,num_kmeans] = size(Cs);
if num_kmeans == 1
    C = Cs;
    return
end

C1 = Cs(:,:,1);
Cs1 = zeros(B,M,num_kmeans);
Cs1(:,:,1) = C1';
for i = 2:num_kmeans
    C = Cs(:,:,i);
    P = permute_endmembers(C1,C);
    C = P*C;
    Cs1(:,:,i) = C';
end

Cs1 = reshape(Cs1, B*M, num_kmeans);
[Cs1, mapping] = pca(Cs1', reduced_dim);

% [log_g, width] = kde(Cs1, Cs1);
for i = 1:size(Cs1,2)
    log_g(:,i) = ksdensity(Cs1(:,i),Cs1(:,i));
end
log_g = prod(log_g,2);

[~,mind] = max(log_g);
C = Cs(:,:,mind);

