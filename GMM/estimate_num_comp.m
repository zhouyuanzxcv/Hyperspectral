function [K,w_jk,mu_jk,sigma_jk,A1] = estimate_num_comp(Y, A, ...
    sizes, shrink_size, max_num_comp, options)
%ESTIMATE_NUM_COMP Summary of this function goes here
%   Detailed explanation goes here
if nargin < 6
    options = [];
end

thresh = 0.99;
ignore_small_prior = 0.1;

[N,M] = size(A);
K = ones(1,M);
w_jk = cell(1,M);
mu_jk = cell(1,M);
sigma_jk = cell(1,M);

A1 = trim_abundance(A, sizes, shrink_size, thresh);

for j = 1:M
    X = Y(A1(:,j), :);
    num = CVIC(X, max_num_comp, options);
    K(j) = num;

    while true
        try % fit_GMM has a small probability to fail
            obj = fit_GMM(X, num, []);
            break;
        catch me
            disp('Failed to perform fit_GMM.');
        end
    end
    w_jk{j} = obj.PComponents;
    mu_jk{j} = obj.mu;
    sigma_jk{j} = obj.Sigma;
end

disp(['The number of components after CVIC are ',num2str(K)]);

%% remove the components that have too smaller priors
for j = 1:M
    X = Y(A1(:,j), :);
    ignore_masks = w_jk{j} < ignore_small_prior;
    w_jk{j}(ignore_masks) = [];
    mu_jk{j}(ignore_masks,:) = [];
    sigma_jk{j}(:,:,ignore_masks) = [];
    
    K(j) = length(w_jk{j});
    obj = struct('PComponents',w_jk{j},'mu',mu_jk{j},'Sigma',sigma_jk{j});
    obj = fit_GMM(X, K(j), obj);
    w_jk{j} = obj.PComponents;
    mu_jk{j} = obj.mu;
    sigma_jk{j} = obj.Sigma;
end

disp(['The number of components after removing small priors are ',num2str(K)]);
disp(['Total number of combinations is ',num2str(prod(K))]);


function num = CVIC(X, max_num_comp, options)
[N,B] = size(X);
if N/2 < B % number of samples is too small to allow 2 components
    num = 1;
    return;
end

warning('off','stats:gmdistribution:FailedToConverge');

% p = randperm(N);
% X = X(p,:);
likelihoods = zeros(N,1);
total_lh = zeros(max_num_comp,1);

num_fold = 5;
% set_size = round(N / num_fold);

for K = 1:max_num_comp
    if (N*(num_fold-1)/num_fold) / K < B % number of samples is too small
        total_lh(K) = -Inf;
        continue;
    end
    for i = 1:num_fold
        X1 = X;
%         indices = (i-1)*set_size+1 : min(i*set_size,N);
        indices = i:num_fold:N;
        test_X = X1(indices,:);
        X1(indices,:) = [];
        try
            obj = fit_GMM(X1, K, []);
            value = calc_log_gmm(test_X, obj.PComponents, obj.mu, obj.Sigma);
            likelihoods(indices) = value;
        catch me
            likelihoods(indices) = -Inf;            
%             disp(['Error in CVIC: ',me.message]);
        end
    end
    total_lh(K) = sum(likelihoods);
end

[max_lh,num] = max(total_lh);
thresh_CVIC = parse_param(options,'thresh_CVIC',[]);
if ~isempty(thresh_CVIC)
    num = find(abs(total_lh - max_lh) <= max_lh * thresh_CVIC, 1);
end

warning('on','stats:gmdistribution:FailedToConverge');

function obj = fit_GMM(X, num, s)
if num > 1
    if isempty(s)
%         s = kmeans(X, num, 'start', 'cluster');
        s = kmeans(X, num);
    end
    options = statset('Display','off');
    obj = gmdistribution.fit(X, num, 'Start', s, 'Options', options);
else % one Gaussian
    [N,B] = size(X);
    obj = [];
    obj.PComponents = 1;
    obj.mu = mean(X,1);
    X1 = X - repmat(obj.mu, N, 1);
    obj.Sigma = (1/N)*(X1'*X1); % same as in EM (not 1/(N-1))
    if N < B % if the number of samples is too small
        obj.Sigma = obj.Sigma + 1e-9*eye(B);
    end
end
    

function A_new = trim_abundance(A, sizes, shrink_size, thresh)
if shrink_size == 0
    A_new = A > thresh;
    return
end

[N,M] = size(A);
A_new = false(N,M);
for j = 1:M
    A1 = reshape(A(:,j)>thresh, sizes);
    se = strel('disk', shrink_size);        
    erodedBW = imerode(A1, se);
    
    shrink_size_r = shrink_size;
    while all(erodedBW(:) == 0) % reduce size of structure element
        shrink_size_r = shrink_size_r - 1;
        se_r = strel('disk', shrink_size_r);
        erodedBW = imerode(A1, se_r);
    end
        
    A_new(:,j) = erodedBW(:);
end

function A_new = trim_abundance1(A, sizes, pure_perc, thresh)
if pure_perc == 1
    A_new = A > thresh;
    return
end

[N,M] = size(A);
A_new = false(N,M);
disk_sizes = (1:100);
for j = 1:M
    A1 = reshape(A(:,j)>thresh, sizes);
    for k = 1:length(disk_sizes)
        se = strel('disk', disk_sizes(k));        
        erodedBW = imerode(A1, se);
        if length(find(erodedBW)) / length(find(A1)) <= pure_perc
            A_new(:,j) = erodedBW(:);
            break;
        end
    end
end
