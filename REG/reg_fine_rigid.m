function options = reg_fine_rigid(rgb1, I1, wl, options)
%REG_FINE_RIGID Summary of this function goes here
%   Detailed explanation goes here

t_start = tic;
% optimize the objective function by block coordinate descent
options = reg_fine_rigid_internal(rgb1, I1, wl, options);

disp(['Elapsed time for fine-scale rigid registration is ',num2str(toc(t_start))]);

function options = reg_fine_rigid_internal(rgb1, I1, wl, options)
I = double(rgb1) / 255;

last_pos = NaN(1,6);

for iter = 1:20
    try
        t_start_iter = tic;
        
        % optimize w.r.t. degree, t, s
        [degree, t, s] = optimize_wrt_T_s(I, I1, wl, options);
        options.degree = degree;
        options.t = t;
        options.s = s;
        
        % optimize w.r.t. lambda
        sigma = optimize_wrt_sigma(I, I1, options);
        options.sigma = sigma;
        
        % check convergence
        if all([degree,t,s,sigma] == last_pos)
            disp(['Iteration ends at ',num2str(iter)]);
            break;
        else
            last_pos = [degree,t,s,sigma];
        end
        disp(['Iteration ',num2str(iter),' lasts ',num2str(toc(t_start_iter))]);
    catch err
        if strcmp(err.identifier, 'Registration:DerivativeBoundary')
            disp(['Stop optimization since boundary touch when ', ...
                'calculating derivatives']);
            break;
        else
            throw(err);
        end
    end
end

    
