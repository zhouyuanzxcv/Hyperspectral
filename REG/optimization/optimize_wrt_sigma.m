function sigma = optimize_wrt_sigma(I, I1, options)
% set the left hand side (warp the hyperspectral image)
options.isRigid = 0;
I1 = transform(I1,eye(3),[1,1],1e-3,size(I1,2),size(I1,1),options);
options.isRigid = 1;

% configure the right hand side (transform the color image)
sigma_step = mean(options.s);
step_size = sigma_step;

pa = ParameterAnalysis();
pa.StepSizeMode = 2;
pa.NeighborhoodMode = 3;
pa.Algorithm = 'BrutalForce';
pa.BrutalForceDepth = 5;
list_params = {'sigma'};
range = [1e-6,4*mean(options.s)]';
fcn_run = @(options1) calc_val(options1,I,I1,options);

options1 = [];
options1.sigma = options.sigma;

options1 = pa.autoParamSelection(fcn_run, options1, list_params, ...
    step_size, range);
sigma = options1.sigma;

function val = calc_val(options1,I,I1,options)
sigma = options1.sigma;
T = create_T(options.degree, options.t);
val = eval_obj_fun(I, I1, T, options.s, sigma, options);
