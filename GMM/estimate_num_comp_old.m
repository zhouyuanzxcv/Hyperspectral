function [K,w_jk,mu_jk,sigma_jk] = estimate_num_comp1(Y, A, w_jk, mu_jk, sigma_jk)
%ESTIMATE_NUM_COMP1 Summary of this function goes here
%   Detailed explanation goes here
max_num_comp = 5;
max_num_configs = 300;

[N,M] = size(A);
K = ones(1,M);
lls = zeros(M, max_num_comp);

for j = 1:M
    X = Y(A(:,j)==1, :);
    lls(j,:) = loglikelihood(X, max_num_comp);
end

lls1 = sort(lls(:),1,'ascend');
lls1(isnan(lls1)) = [];
for i = 1:length(lls1)
    thresh = lls1(i);
    for j = 1:M
        num = find(lls(j,:)>=thresh, 1, 'first');
        if isempty(num)
            num = max_num_comp;
        end
        K(j) = num;
    end
    if prod(K) > max_num_configs
        thresh = lls1(i-1);
        break;
    end
end

% thresh = min(lls(:,end));
for j = 1:M
    X = Y(A(:,j)==1, :);

    num = find(lls(j,:)>=thresh, 1, 'first');
    if isempty(num)
        num = max_num_comp;
    end
    
    K(j) = num;
    
    s = kmeans(X,num);
    options = statset('Display','off');
    obj = gmdistribution.fit(X, num, 'Start', s, 'Options', options);
    w_jk{j} = obj.PComponents;
    mu_jk{j} = obj.mu;
    sigma_jk{j} = obj.Sigma;
    
    disp(['The number of components for endmember ',num2str(j), ...
        ' is estimated to be ',num2str(num)]);
end

function lls = loglikelihood(X, max_num_comp)
lls = zeros(1, max_num_comp);
[N,B] = size(X);
% [log_g, s] = kde(X, X);

for k = 1:max_num_comp
    try
        % use kmean to initialize the GMM is important
        s = kmeans(X,k);
        options = statset('Display','off');
        obj = gmdistribution.fit(X, k, 'Start', s, 'Options', options);
        lls(k) = (-1/N) * obj.NlogL;
    catch me
        lls(k) = NaN;            
        disp(['Error in loglikelihood: ',me.message]);
    end
end

function num = BIC(X)
max_num_comp = 5;

divs = zeros(1, max_num_comp);
BICs = zeros(1, max_num_comp);
AICs = zeros(1, max_num_comp);

[log_g, s] = kde(X, X);

for k = 1:max_num_comp
    % use kmean to initialize the GMM is important
    s = kmeans(X,k);
    options = statset('Display','final');
    obj = gmdistribution.fit(X, k, 'Start', s, 'Options', options);

    divs(k) = KL_div(log_g, X, obj.PComponents, obj.mu, obj.Sigma);
    BICs(k) = obj.BIC;
    AICs(k) = obj.AIC;
end

[~,num] = min(divs);
[~,num] = min(BICs);
[~,num] = min(AICs);

function num = empirical_num_comp(X)
[mappedX, mapping] = pca(X,2);
X1 = mappedX(:,1);
n = histc(X1,linspace(min(X1),max(X1),30));
h = fspecial('gaussian',2*6+1,2);
n1 = imfilter(n,h);
ind = find(n1 > [n1(2:end);n1(end)] & n1 > [n1(1);n1(1:end-1)]);
num = length(ind);