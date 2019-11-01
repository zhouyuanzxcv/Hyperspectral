function [I1,wl1,rgb1,g,H1,extra,s,G,S] = estimate_PSF_SRF(I,wl,bbl,rgb,s,extra,options)
%ESTIMATE_PSF_SRF Estimate point spread function (PSF) and spectral 
%response function (SRF)
disp('Start estimating PSF and SRF');
show_fig = parse_param(options,'show_fig',0);

I1 = I(:,:,bbl);
wl1 = wl(bbl);
rgb1 = rgb;

Y = reshape_hsi(I1);
X = reshape(rgb1, [size(rgb1,1)*size(rgb1,2), size(rgb1,3)]);
[N,B] = size(Y);
b = size(X,2);
R2 = s + extra*2;
R = prod(R2);

convergence_t = 1e-3;

SRF_range = parse_param(options,'SRF_range','color');
% construct S
[I_sel,bands_sel] = select_relevant_bands(I1,wl1,SRF_range);
wl_sel = wl1(bands_sel)';
S = zeros(B,length(wl_sel));
for i = 1:size(S,2)
    S(bands_sel(i),i) = 1;
end

% construct C
Cs = construct_C(I1,rgb1,s,extra);
C = cat(2,Cs{:});

Y1 = [ones(N,1),Y*S];

% construct D
D = zeros(N*b,R);
for i = 1:N
    D((i-1)*b+1:i*b,:) = X'*Cs{i}';
end

% contruct options
options = [];
options.D = D;
options.Y1 = Y1;
options.D1D = D'*D;
options.D1YI = D'*kron(Y1,eye(b));

g = (1/R) * ones(R,1);
G = g2G(C,g,N);


delta_t0 = 1e-4;
delta_t_g = delta_t0;

errors = [];
for iter = 1:200
    % Update H1
    H1 = solve_for_H1(G*X, Y1, struct('lambda',1e-4));
    options.H1 = H1;
    
    % Update G
    der_g = calc_der_g(g, options);
    options.der_g = der_g;
    delta_t_g = calc_time_step_adaptive(@eval_obj_fun_g, @update_g, ...
        g, options, delta_t_g, delta_t0);
    g = update_g(g, options, delta_t_g);
    options.g = g;
    G = g2G(C,g,N);
    
    % Test convergence
    errors(end+1) = eval_obj_fun(options);

    if test_convergence(errors, convergence_t)
        break;
    end
    
    if mod(iter,100) == 0
        disp(['Process iteration ',num2str(iter)]);
    end
end

if show_fig
    figure('name','PSF'), mesh(reshape(g,R2));
    figure('name','SRF'), plot(H1);
end

function der_g = calc_der_g(g, options)
vecH1 = (options.H1)';
vecH1 = vecH1(:);
der_g = options.D1D * g - options.D1YI * vecH1;

function val = eval_obj_fun_g(g, options)
vecYH = options.Y1 * options.H1;
vecYH = vecYH';
vecYH = vecYH(:);
val = sum((options.D * g - vecYH).^2);

function g_new = update_g(g, options, delta_t)
der_g = options.der_g;

g_new = g - delta_t * der_g;
g_new = project_to_simplex(g_new');
g_new = g_new';

function val = eval_obj_fun(options)
g = options.g;
val = eval_obj_fun_g(g,options);

