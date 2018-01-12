function [mu_jk,sigma_jk,w_jk,K,A] = gmm_init(I1, M, options)
%INIT_GMM Summary of this function goes here
%   Detailed explanation goes here
sigma0 = 1e-2;

if isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end

[Y,~,rows,cols] = reshape_hsi(I1,[]);
[N,B] = size(Y);

C = find_init_clusters_by_kmeans(I1,M);
% [idx,C] = kmeans(Y,M,'start',C);

mu_jk = cell(1,M);
sigma_jk = cell(1,M);
w_jk = cell(1,M);

K = ones(1,M);

for j = 1:M
    mu_jk{j}(1,:) = C(j,:);
%     X1 = Y(idx==j,:);
%     X1_c = mean(X1,1);
%         X1_rc = X1 - repmat(X1_c, size(X1,1), 1);
%         sigma_jk{j}(:,:,1) = (X1_rc' * X1_rc)/size(X1,1);
    sigma_jk{j}(:,:,1) = sigma0^2 * eye(B);
    w_jk{j}(1,1) = 1;
end

R = C;
A = Y*R'*inv(R*R'+eye(M)*1e-6);
A = project_to_simplex(A);

%% show the initial condition
names = cell(1,M);
for i = 1:M
    names{i} = ['endmember ',num2str(i)];
end

% show scatter plot
endmember_scatter_plot_end_var(Y,w_jk,mu_jk,sigma_jk,names,options);
set(gcf,'name','Scatter plot of the estimated initial Gaussians');

show_abundances(A,rows,cols);
set(gcf,'name','Initial abundances');

function C = find_init_clusters_by_kmeans(I1,M)
disp('Start kmeans to find the initial clusters ...');

if 0
    gauss_sigma = 1.5;
    h = fspecial('gaussian', 2*round(3*gauss_sigma)+1, gauss_sigma);
    for i = 1:size(I1,3)
        I2 = I1(:,:,i);
        I1(:,:,i) = imfilter(I2,h,'replicate');
    end
end

[Y,~,rows,cols] = reshape_hsi(I1,[]);
[N,B] = size(Y);

if 1
    tic
    num_kmeans1 = 1;
    Css = zeros(B,M,num_kmeans1);
    for j = 1:num_kmeans1
        num_kmeans = 200;
        Cs = zeros(B,M,num_kmeans);
        for i = 1:num_kmeans
            Y1 = Y;
    %         Y1 = Y(rand(size(Y,1),1) < 0.1, :);
            [~, C] = kmeans(Y1, M);
            Cs(:,:,i) = C';
        end

        Css(:,:,j) = find_most_frequent(Cs, 5);
    end
    C = find_most_frequent(Css, 2);
    C = C';
    toc
else
    tic
    [idx,C] = kmeans(Y,M,'Replicates',100);
    toc
end

function C = find_most_frequent(Cs, reduced_dim)
% reduced_dim = 5;
[B,M,num_kmeans] = size(Cs);
if num_kmeans == 1
    C = Cs;
    return
end

C1 = Cs(:,:,1);
for i = 2:num_kmeans
    C = Cs(:,:,i);
    P = permute_endmembers(C1',C');
    C = C*P';
    Cs(:,:,i) = C;
end

Cs = reshape(Cs, B*M, num_kmeans);
if 1
    [Cs1, mapping] = pca(Cs', reduced_dim);
else
    nn = max(num_kmeans/20, 2);
    Cs1 = LE(Cs','nn',nn,Inf,reduced_dim);
end
[log_g, width] = kde(Cs1, Cs1);
[~,mind] = max(log_g);
C = Cs(:,mind);
C = reshape(C,B,M);


function [A, mu_jk, sigma_jk, w_jk, K] = segment_and_estimate_from_A(I1, M)
options = [];

[Y,~,rows,cols] = reshape_hsi(I1,[]);
[N,B] = size(Y);

% Use GMM on pure pixel to estimate K and initialize mu and sigma
[A, thresh] = estimate_A(I1, M);
[mu_jk, sigma_jk, w_jk, K, Xs] = initialize_from_A(Y, A, thresh, M);
options.Xs = Xs;

% show histogram
hist_end_var(Y,A,names,thresh);
set(gcf,'name','Histogram of picked pure pixels for initial K estimation');


function [mu_jk, sigma_jk, w_jk, K, Xs] = initialize_from_A(Y, A, thresh, M)
mu_jk = cell(1,M);
sigma_jk = cell(1,M);
w_jk = cell(1,M);

K = ones(1,M);
Xs = cell(1,M);
for j = 1:M
    ind = find(A(:,j) >= max(A(:,j))*thresh);
    X = Y(ind,:);
    Xs{j} = X;
    if 1 % use GMM in matlab on the pure pixels
        if 0
            max_components = 4;
            options = statset('Display','final');
            AIC = zeros(1,max_components);
            BIC = zeros(1,max_components);
            objs = cell(1,max_components);
            for k = 1:max_components
                try
                    s = kmeans(X,k);
                    objs{k} = gmdistribution.fit(X, k, 'Start', s, 'Options', options);
                    AIC(k) = objs{k}.AIC;
                    BIC(k) = objs{k}.BIC;
                catch me
                    disp(me.message);
                    objs{k} = [];
                    AIC(k) = Inf;
                    BIC(k) = Inf;
                end
            end
            [~,min_comp] = min(AIC);

            K(j) = min_comp;
            obj = objs{min_comp};
            w_jk{j} = obj.PComponents;
            mu_jk{j} = obj.mu;
            sigma_jk{j} = obj.Sigma;
        else
            k = empirical_num_comp(X);
            K(j) = k;
            % use kmean to initialize the GMM is important
            s = kmeans(X,k);
            options = statset('Display','final');
            obj = gmdistribution.fit(X, k, 'Start', s, 'Options', options);
            w_jk{j} = obj.PComponents;
            mu_jk{j} = obj.mu;
            sigma_jk{j} = obj.Sigma;
        end
    else % use k-means in matlab on the pure pixels
        [idx,C] = kmeans(X,K(j));
        mu_jk{j} = C;
        for k = 1:K(j)
            X1 = X(idx==k,:);
            X1_c = mean(X1,1);
            X1_rc = X1 - repmat(X1_c, size(X1,1), 1);
%             sigma_jk{j}(:,:,k) = (X1_rc' * X1_rc)/size(X1,1);
            sigma_jk{j}(:,:,k) = 0.5 * (X1_rc' * X1_rc)/size(X1,1);
            w_jk{j}(1,k) = size(X1,1) / size(X,1);
        end
    end
end


function [A, thresh] = estimate_A(I1, M)
[Y,~,rows,cols] = reshape_hsi(I1,[]);
[N,B] = size(Y);

if 0
    options.beta1 = 0.01;
    options.beta2 = 0;
    options.rho1 = 0.05; % 0.005, 0.1
    [A,R,mu,sigma,var_dirs,var_amts] = scm(I1,M,options);

    show_abundances(A,rows,cols);
    
    thresh = 0.95;

    % [idx,C] = kmeans(Y,M);
    % A = zeros(N,M);
    % for j = 1:M
    %     A(idx==j,j) = 1;
    % end
else
    % Use graph cut to segment the image
    A = zeros(N,M);
    
    [L,C] = PegahGraphCut(I1,M);
    colors = distinguishable_colors(M);
    figure('name','Segmentation results');
    imshow(uint8(L),colors);
    cbh = colorbar;
    keep_perc = 0.7;
    for j = 1:M
        binary = L==j-1;
        total = length(find(binary));
        for s = 1:10
            se = strel('disk',s);
            eroded = imerode(binary,se);
            if length(find(eroded)) <= total * keep_perc
                break;
            end
        end
        A(eroded(:)==1,j) = 1;
    end
    thresh = 1;
    
    R = zeros(M,B);
    for j = 1:M
        R(j,:) = mean(Y(A(:,j)==1,:),1);
    end
    A1 = Y*R'*inv(R*R'+eye(M)*1e-6);
    A1 = project_to_simplex(A1);
    for j = 1:M
        idx = A(:,j)==1;
        A1(idx,:) = A(idx,:);
    end
    A = A1;
end

