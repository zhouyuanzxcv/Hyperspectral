function options = reg_fine_nonrigid(rgb1, I1, wl, options)
%REG_FINE_NONRIGID Summary of this function goes here
%   Detailed explanation goes here
t_start = tic;

I = double(rgb1) / 255;

nonrigid_model = parse_param(options,'nonrigid_model','LSQ freeform');

% used for B-spline with MI metric
O_trans = [];
Spacing = [];

% old variables to check convergence
last_pos = NaN(1,6);
U_old = options.U;
V_old = options.V;

fix_rigid_params = 0;
if isfield(options,'fix_rigid_params') && options.fix_rigid_params
    fix_rigid_params = 1;
end

for iter = 1:20
    t_start_iter = tic;
    % optimize w.r.t. degree, t, s
    if fix_rigid_params
        disp('Skip rigid transformation of color image');
        degree = options.degree;
        t = options.t(:)';
        s = options.s(:)';
    else
        [degree, t, s] = optimize_wrt_T_s(I, I1, wl, options);
        options.degree = degree;
        options.t = t;
        options.s = s;
    end

    % optimize w.r.t U,V
    if strcmp(nonrigid_model,'LSQ freeform')
        [U,V] = optimize_wrt_V(rgb1, I1, wl, options);
    elseif strcmp(nonrigid_model,'MI bspline')
        [U,V,O_trans,Spacing] = optimize_wrt_V_by_MI_spline(rgb1, I1, wl, ...
            O_trans, Spacing, options);
    else
        disp('No nonrigid model specified');
    end
    options.U = U;
    options.V = V;

    % optimize w.r.t. lambda
    sigma = optimize_wrt_sigma(I, I1, options);
    options.sigma = sigma;
    
    % test convergence. When max is less than 0.1, mean is less than 0.01
    if all([degree,t,s,sigma] == last_pos) && max(abs(U(:) - U_old(:))) ...
            < 0.1 && max(abs(V(:) - V_old(:))) < 0.1
        disp(['Iteration ends at ',num2str(iter)]);
        break;
    else
        last_pos = [degree,t,s,sigma];
        U_old = U;
        V_old = V;
    end
    disp(['Iteration ',num2str(iter),' lasts ',num2str(toc(t_start_iter))]);
    
    % fix rigid parameters after some iterations to avoid image drifting
    if iter >= 1
        fix_rigid_params = 1;
    end
end

disp(['Elapsed time for nonrigid registration is ',num2str(toc(t_start))]);


function [U,V,O_trans,Spacing] = optimize_wrt_V_by_MI_spline(rgb1, I1, ...
    wl, O_trans, Spacing, options)
degree = options.degree;
t = options.t;
s = options.s;
T = create_T(degree, t);
sigma = options.sigma;

Istatic = transform(double(rgb1(:,:,1)),T,s,sigma,size(I1,2),size(I1,1),options);

rgb_ind = find_rgb_ind(wl);
Imoving= I1(:,:,rgb_ind(1)) * 255;

reg_opts = [];
reg_opts.Similarity = 'mi';
reg_opts.Registration = 'NonRigid';

% Try to use the last time's B-spline parameters as input for the
% current iteration.
if ~isempty(O_trans) && ~isempty(Spacing) % not the first time run
    reg_opts.Grid = O_trans;
    reg_opts.Spacing = Spacing; % spacing is already the minimum value 2
    reg_opts.MaxRef = 0;
end

[Ireg,O_trans,Spacing,M,B,F] = image_registration(Imoving,Istatic,reg_opts);
% show_nonrigid_result(Imoving,[],B(:,:,2),B(:,:,1));
U = B(:,:,2);
V = B(:,:,1);
