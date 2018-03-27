function [T,degree,t,s,sigma,U,V] = reg_hyper_rgb(rgb1, I1, wl, options)
%REG_HYPER_RGB Summary of this function goes here
%   Detailed explanation goes here
t_start_all = tic;

% set default parameters
options.lambda = parse_param(options,'lambda',1e-3);

%% prepare calculating objective function
% options.L = calc_linear_operator(I1, wl);
% K = 3;
% [J1,J2] = find_nearest_farthest_points(reshape_hsi(I1), K);
% options.J1 = J1;
% options.J2 = J2;
% [L1,L2] = calc_difference_operator(reshape_hsi(I1), K);
% options.L1 = L1;
% options.L2 = L2;

options = prepare_calc_reg_obj(I1,wl,options);

options.U = zeros(size(I1,1), size(I1,2));
options.V = options.U;

options.isRigid = 1; % used in transform
options.useTranslationField = 1;

options.wl = wl;

options.show_figure = parse_param(options,'show_figure',0);

%% find initial registration by phase correlation
if isfield(options,'fix_rigid_params') && options.fix_rigid_params
    disp('Skip initialization of rigid parameters');
    options.sigma = 2 * mean(options.s);
else
    if 1
        [degree,t,s,sigma] = reg_init(rgb1, I1, wl, options);
        save('result_init.mat','degree','t','s','sigma');
    else
        load('result_init.mat');
    end

    options.degree = degree;
    options.s = s;
    options.t = t;
    options.sigma = sigma;
    disp('Initial registration gives ');
    options
end

if options.show_figure % show initial condition in the fine scale
    I4 = transform(double(rgb1),create_T(degree,t),s,sigma,size(I1,2),...
        size(I1,1),options);
    figure('name','initial condition in fine scale');
    imshow(uint8(I4));
end

%% find rigid/nonrigid transformation in fine scale
reg_method = parse_param(options, 'reg_method', 'nonrigid');
if strcmp(reg_method, 'rigid')
    % rigid case
    options = reg_fine_rigid(rgb1, I1, wl, options);    
elseif strcmp(reg_method, 'nonrigid')
    % nonrigid case
    options = reg_fine_nonrigid(rgb1, I1, wl, options);    
end

disp('Fine-scale registration gives');
options

degree = options.degree;
t = options.t;
s = options.s;
T = create_T(options.degree, options.t);
sigma = options.sigma;
U = options.U;
V = options.V;
disp(['The total registration time is ',num2str(toc(t_start_all))]);


% function U = calc_U(I, wl)
% Y1 = calc_Y1(I, wl);
% [U,S,V] = svd(Y1,0);
% % U1 = U(:,size(Y1,2)+1:end);
% % U1 = U1';
% 

function [J1,J2] = find_nearest_farthest_points(X, K)
step = 100;
n = size(X,1);
J1 = cell(n,1);
J2 = cell(n,1);

for i1 = 1:step:n
    i2 = i1 + step - 1;
    if (i2 > n)
        i2 = n;
    end;

    XX = X(i1:i2,:);
    D = calc_distance_matrix(XX,X);
    [Z,I] = sort(D,2);

    Z = Z(:,2:end);
    I = I(:,2:end);
    for i = i1:i2
        J1{i} = I(i-i1+1,1:K);
        J2{i} = I(i-i1+1,size(Z,2)-K+1:end);
    end
end

function [L1,L2] = calc_difference_operator(Y, K)
[N,B] = size(Y);
[J1,J2] = find_nearest_farthest_points(Y, K);
J1 = cell2mat(J1);
J2 = cell2mat(J2);
L1 = cell(1,K);
L2 = cell(1,K);

for k = 1:K
    J_1 = [J1(:,k);(1:N)'];
    J_2 = [J2(:,k);(1:N)'];
    I = [(1:N)';(1:N)'];
    V = [ones(N,1);-ones(N,1)];
    L1{k} = sparse(I,J_1,V,N,N);
    L2{k} = sparse(I,J_2,V,N,N);
end

% function L = calc_linear_operator(I, wl)
% % calculate neighborhood based on 8-neighbor system
% epsilon = 0.01;
% K = 8;
% [rows,cols,B] = size(I);
% Y_ori = reshape_hsi(I);
% rgb_ind = find_rgb_ind(wl);
% srf_width = 2;
% Y = [];
% for ind = rgb_ind
%     Y = [Y, Y_ori(:, ind - srf_width : ind + srf_width)];
% end
% 
% N = rows * cols;
% N1 = (rows-2) * (cols-2);
% J2 = cell(1,N1);
% A2 = cell(1,N1);
% ind_J2 = 1;
% for j = 2:cols-1
%     for i = 2:rows-1
%         y = [i-1,i,i+1,i-1,i+1,i-1,i,i+1];
%         x = [j-1,j-1,j-1,j,j,j+1,j+1,j+1];
%         inds = sub2ind([rows,cols],y,x);
%         curr = sub2ind([rows,cols],i,j);
%         Z = repmat(Y(curr,:),[K,1]) - Y(inds,:);
%         C = Z*Z' + epsilon*eye(K);
%         wn = C\ones(K,1);
%         wn = wn/sum(wn);
%         J2{ind_J2} = [inds,curr]';
%         A2{ind_J2} = [wn;-1];
%         ind_J2 = ind_J2 + 1;
%     end
% end
% 
% L = cellarr2sparse(J2,A2,N1,N);
