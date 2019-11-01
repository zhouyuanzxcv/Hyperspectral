function R = optimize_gradient_descent(I, I1, options)
delta_t0 = parse_param(options,'delta_t0',1e-6);
max_iter = parse_param(options,'max_iter',500);
convergence_t = parse_param(options,'convergence_t',0.001);

R = init_R(I,I1,options);

errors = eval_obj_fun_R(R, options);

delta_t_R = delta_t0;

for iter = 1:max_iter
    der_R = calc_der_R(R, options);
    options.der_R = der_R;
    delta_t_R = calc_time_step_adaptive(@eval_obj_fun_R, @update_R, ...
        R, options, delta_t_R, delta_t0);
    R = update_R(R, options, delta_t_R);
    
    % test convergence
    errors(end+1) = eval_obj_fun_R(R, options);
    if test_convergence(errors, convergence_t)
        break;
    end

end

function R = init_R(I,I1,options)
[rows1,cols1,b] = size(I1);
[rows,cols,B] = size(I);
N = rows*cols;
extra = options.extra;
sy = (rows1 - 2*extra(2)) / rows; % should be integer
sx = (cols1 - 2*extra(1)) / cols;

R = zeros(rows1*cols1, B);
for i = 1:N
    y = options.Y(i,:);
    
    [r,c] = ind2sub([rows,cols],i);
    x1_inds = (c-1)*sx+1:c*sx;
    y1_inds = (r-1)*sy+1:r*sy;
    [X,Y] = meshgrid(x1_inds+extra(1), y1_inds+extra(2));
    inds = sub2ind([rows1,cols1],Y(:),X(:));
    R(inds,:) = repmat(y, [length(inds),1]);
end 


function der_R = calc_der_R(R, options)
beta = options.beta;
gamma = options.gamma;
F = options.F;
G = options.G;
Z = options.Z;
L = options.L;

der_R = gamma*G'*(G*R) + beta*L*R + (1-gamma)*(R*F)*F' - Z;


function R_new = update_R(R, options, delta_t)
der_R = options.der_R;

R_new = R - delta_t * der_R;
R_new(R_new < 0) = 0;


function val = eval_obj_fun_R(R, options)
beta = options.beta;
gamma = options.gamma;
G = options.G;
Y = options.Y;
F = options.F;
X1 = options.X1;
L = options.L;

val = gamma*sum(sum((G*R-Y).^2)) + (1-gamma)*sum(sum((R*F-X1).^2));
val = val + beta*trace(R'*L*R);