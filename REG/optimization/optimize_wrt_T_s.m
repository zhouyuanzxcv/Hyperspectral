function [degree, t, s] = optimize_wrt_T_s(I, I1, wl, options)
% warp the hyperspectral image based on options.U and options.V
options.isRigid = 0;
I1 = transform(I1,eye(3),[1,1],1e-3,size(I1,2),size(I1,1),options);
options.isRigid = 1;

options = prepare_calc_reg_obj(I1,wl,options);

% optimize w.r.t degree, t and s
[degree,t,s] = optimize_by_bcd(I, I1, options);
% [degree,t,s] = optimize_by_brutal_force(I, I1, options);
% [degree,t,s] = optimize_by_gradient_descent(I, I1, options);


%--------------------- objective function evaluation ---------------------%
function val = calc_val(options1,I,I1,options)
degree = options1.degree;
t = [options1.tx, options1.ty];
s = [options1.sx, options1.sy];
T = create_T(degree,t);
sigma = options.sigma;
val = eval_obj_fun(I, I1, T, s, sigma, options);


%--------------------- block coordinate descent --------------------------%
function [degree,t,s] = optimize_by_bcd(I, I1, options)
%% old implementation
% t_step = 0.2 * mean(options.s);
% s_step = 0.02 * mean(options.s);
% step_size = [0.2, t_step, t_step, s_step, s_step];
% 
% options.degree = optimze_wrt_degree(I,I1,options,step_size(1));
% options.t(1) = optimze_wrt_tx(I,I1,options,step_size(2));
% options.t(2) = optimize_wrt_ty(I,I1,options,step_size(3));
% options.s(1) = optimize_wrt_sx(I,I1,options,step_size(4));
% options.s(2) = optimize_wrt_sy(I,I1,options,step_size(5));
% 
% degree = options.degree;
% t = options.t;
% s = options.s;

%% new implementation
t_step = 0.1 * mean(options.s);
s_step = 0.01 * mean(options.s);
step_size = [0.1, t_step, t_step, s_step, s_step];

pa = ParameterAnalysis();
pa.StepSizeMode = 2; % addition
pa.NeighborhoodMode = 3*ones(1,5);

% pa.MaxNumberOfIterations = 1;
pa.Algorithm = 'BlockCoordinateDescent';
pa.BlockCoordinateDescentDepth = 5;
pa.BlockCoordinateDescentGroup = {1,[2,3],[4,5]};

list_params = {'degree','tx','ty','sx','sy'};
range = [-1e6,1e6;-1e6,1e6;-1e6,1e6;-1e6,1e6;-1e6,1e6]';
fcn_run = @(options1) calc_val(options1,I,I1,options);

options1 = [];
options1.degree = options.degree;
options1.tx = options.t(1);
options1.ty = options.t(2);
options1.sx = options.s(1);
options1.sy = options.s(2);

options1 = pa.autoParamSelection(fcn_run, options1, list_params, ...
    step_size, range);
degree = options1.degree;
t = [options1.tx, options1.ty];
s = [options1.sx, options1.sy];


function val = optimize_wrt_single_variable(options,step_size,name,fcn_run)
pa = ParameterAnalysis();
pa.StepSizeMode = 2;
pa.NeighborhoodMode = 3;
pa.Algorithm = 'BrutalForce';
pa.BrutalForceDepth = 5;
list_params = {name};
range = [-1e6,1e6]';

options1 = [];
switch name
    case 'degree'
        options1.(name) = options.(name);
    case 'tx'
        options1.(name) = options.t(1);
    case 'ty'
        options1.(name) = options.t(2);
    case 'sx'
        options1.(name) = options.s(1);
    case 'sy'
        options1.(name) = options.s(2);
    otherwise
end

options1 = pa.autoParamSelection(fcn_run, options1, list_params, ...
    step_size, range);
val = options1.(name);


function degree = optimze_wrt_degree(I,I1,options,step_size)
fcn_run = @(options1) calc_val_degree(options1,I,I1,options);
degree = optimize_wrt_single_variable(options,step_size,'degree',fcn_run);

function tx = optimze_wrt_tx(I,I1,options,step_size)
fcn_run = @(options1) calc_val_tx(options1,I,I1,options);
tx = optimize_wrt_single_variable(options,step_size,'tx',fcn_run);

function ty = optimize_wrt_ty(I,I1,options,step_size)
fcn_run = @(options1) calc_val_ty(options1,I,I1,options);
ty = optimize_wrt_single_variable(options,step_size,'ty',fcn_run);

function sx = optimize_wrt_sx(I,I1,options,step_size)
fcn_run = @(options1) calc_val_sx(options1,I,I1,options);
sx = optimize_wrt_single_variable(options,step_size,'sx',fcn_run);

function sy = optimize_wrt_sy(I,I1,options,step_size)
fcn_run = @(options1) calc_val_sy(options1,I,I1,options);
sy = optimize_wrt_single_variable(options,step_size,'sy',fcn_run);

function val = calc_val_degree(options1,I,I1,options)
degree = options1.degree;
T = create_T(degree, options.t);
sigma = options.sigma;
val = eval_obj_fun(I, I1, T, options.s, sigma, options);

function val = calc_val_tx(options1,I,I1,options)
tx = options1.tx;
t = options.t;
t(1) = tx;
T = create_T(options.degree, t);
sigma = options.sigma;
val = eval_obj_fun(I, I1, T, options.s, sigma, options);

function val = calc_val_ty(options1,I,I1,options)
ty = options1.ty;
t = options.t;
t(2) = ty;
T = create_T(options.degree, t);
sigma = options.sigma;
val = eval_obj_fun(I, I1, T, options.s, sigma, options);

function val = calc_val_sx(options1,I,I1,options)
sx = options1.sx;
s = options.s;
s(1) = sx;
T = create_T(options.degree, options.t);
sigma = options.sigma;
val = eval_obj_fun(I, I1, T, s, sigma, options);

function val = calc_val_sy(options1,I,I1,options)
sy = options1.sy;
s = options.s;
s(2) = sy;
T = create_T(options.degree, options.t);
sigma = options.sigma;
val = eval_obj_fun(I, I1, T, s, sigma, options);


%--------------------- brutal force search -------------------------------%
function [degree,t,s] = optimize_by_brutal_force(I, I1, options)
t_step = 0.2 * mean(options.s);
s_step = 0.02 * mean(options.s);
step_size = [0.2, t_step, t_step, s_step, s_step];

pa = ParameterAnalysis();
pa.StepSizeMode = 2;
pa.NeighborhoodMode = 3*ones(1,5);
pa.Algorithm = 'BrutalForce';
pa.BrutalForceDepth = 5;
list_params = {'degree','tx','ty','sx','sy'};
range = [-1e6,1e6;-1e6,1e6;-1e6,1e6;-1e6,1e6;-1e6,1e6]';
fcn_run = @(options1) calc_val(options1,I,I1,options);

options1 = [];
options1.degree = options.degree;
options1.tx = options.t(1);
options1.ty = options.t(2);
options1.sx = options.s(1);
options1.sy = options.s(2);

options1 = pa.autoParamSelection(fcn_run, options1, list_params, ...
    step_size, range);
degree = options1.degree;
t = [options1.tx, options1.ty];
s = [options1.sx, options1.sy];


%--------------------- gradient descent ---------------------------------%
function [degree,t,s] = optimize_by_gradient_descent(I, I1, options)
delta_t0 = 1e-4;
delta_t = delta_t0;
sigma = options.sigma;
T_s = [options.degree; options.t; options.s];

errors = eval_obj_fun(I, I1, create_T(T_s(1), T_s(2:3)), T_s(4:5), ...
    sigma, options);
for iter = 1:100
    der_T_s = calc_der_T_s(I, I1, T_s, sigma, options);
    options.der_T_s = der_T_s;
    eval_obj_fun_T_s = @(T_s, options) eval_obj_fun(I, I1, ...
        create_T(T_s(1), T_s(2:3)), T_s(4:5), sigma, options);
    delta_t = calc_time_step_adaptive(eval_obj_fun_T_s, @update_T_s, ...
        T_s, options, delta_t, delta_t0);
    T_s = update_T_s(T_s, options, delta_t);

    errors = [errors, eval_obj_fun(I, I1, create_T(T_s(1), ...
        T_s(2:3)), T_s(4:5), sigma, options)];
    
    if test_convergence(errors, 1e-2)
        break;
    end
end
degree = T_s(1);
t = T_s(2:3);
s = T_s(4:5);

function T_s = update_T_s(T_s, options, delta_t)
T_s = T_s - delta_t * options.der_T_s;
T_s(1) = mod(T_s(1), 360);
s = T_s(4:5);
s(s < 0) = 1e-4;
T_s(4:5) = s;

function der_T_s = calc_der_T_s(I, I1, T_s, lambda, options)
der_T_s = zeros(length(T_s), 1);
step = 0.1;

for i = 1:length(T_s)
    prev = T_s;
    forw = T_s;
    prev(i) = T_s(i) - step;
    forw(i) = T_s(i) + step;
    val_prev = eval_obj_fun(I, I1, create_T(prev(1), prev(2:3)), ...
        prev(4:5), lambda, options);
    val_forw = eval_obj_fun(I, I1, create_T(forw(1), forw(2:3)), ...
        forw(4:5), lambda, options);
    if (isinf(val_prev) || isinf(val_forw))
        throw(MException('Registration:DerivativeBoundary', ['Touch ', ...
            'boundary when calculating derivatives.']));
    end
    der_T_s(i) = (val_forw - val_prev) / (2*step);
end
