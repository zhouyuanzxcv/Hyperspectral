function [mu_jk,sigma_jk,w_jk,K,A] = gmm_init(I1, M, options)
%INIT_GMM Summary of this function goes here
%   Detailed explanation goes here
sigma0 = 0.1;

if isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end

[Y,~,rows,cols] = reshape_hsi(I1,[]);
[N,B] = size(Y);

C = find_init_clusters_by_kmeans(I1,M);

mu_jk = cell(1,M);
sigma_jk = cell(1,M);
w_jk = cell(1,M);

K = ones(1,M);

for j = 1:M
    mu_jk{j}(1,:) = C(j,:);
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

[Y,~,rows,cols] = reshape_hsi(I1,[]);

tic
num_kmeans = 100;
C = kmeans_stable(Y,M,num_kmeans);
toc

