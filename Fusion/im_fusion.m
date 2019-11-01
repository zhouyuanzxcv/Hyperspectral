function I2 = im_fusion(I,wl,I1,options)
%IM_FUSION Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   I - rows by cols by B hyperspectral image, in the range [0 1]
%   wl - 1 by B wavelength
%   I1 - rows1 by cols1 by b RGB image, in the range [0 255] or [0 1]
if nargin < 4
    options = [];
end

t_start = tic;

I1 = double(I1);
scale_rgb = 0;
if max(I1(:)) > 10 % in the range 0 - 255
    I1 = I1/255;
    scale_rgb = 1;
end

has_reg_params = parse_param(options,'has_reg_params',0);
bbl = parse_param(options,'bbl',true(1,size(I,3)));
show_fig = parse_param(options,'show_fig',0);
has_PSF_SRF = parse_param(options,'has_PSF_SRF',0);
regularization = parse_param(options,'regularization','LLMM'); %LMM or LLMM
sylvester_solver_mode = parse_param(options,'sylvester_solver_mode','adaptive_accurate');

%% automatically infer s and extra and check size match
[rows,cols,B] = size(I);
[rows1,cols1,b] = size(I1);

s = parse_param(options,'s', [floor(cols1/cols),floor(rows1/rows)]);
extra = parse_param(options,'extra',max(([cols1,rows1]-s.*[cols,rows])/2, [1,1]));

pad_color = 0;
if all([cols,rows].*s + extra*2 == [cols1,rows1])
    disp('Image sizes match exactly');
elseif all([cols,rows].*s == [cols1,rows1])
    % add padding for extra margin
    I1 = pad_rgb_margin(I1, extra);
    disp('Pad color image with extra margin');
    [rows1,cols1,b] = size(I1);
    pad_color = 1;
else
    error('Image sizes do not match');
end

%% decode PSF and SRF
if has_PSF_SRF
    H1 = options.H1;
    G = options.G;
    S = options.S;
    s = options.s;
    g = options.g;
    extra = options.extra;
    if scale_rgb
        H1 = H1/255;
    end
    I(:,:,~bbl) = [];
    wl(~bbl) = [];
elseif has_reg_params
    reg = options.reg_params;
    % bad bands are removed in interpret_reg_params
    [I,wl,I1,g,H1,extra,s,G,S] = interpret_reg_params(I,wl,bbl,...
        I1,reg.U2,reg.V2,reg.rho,reg.s2,reg.T2,reg.sigma2);
else    
    % bad bands are removed in estimate_PSF_SRF
    [I,wl,I1,g,H1,extra,s,G,S] = estimate_PSF_SRF(I,wl,bbl,I1,s,extra,options);
end

%% adjust parameters
[rows,cols,B] = size(I);
[rows1,cols1,b] = size(I1);

N = rows * cols;
N1 = rows1 * cols1;

h0 = H1(1,:)';
H = H1(2:end,:);
Y = reshape_hsi(I);
X = reshape(I1, [N1, b]);

beta = parse_param(options,'beta',1);
% beta = max(beta, 1e-12);
% 
% beta2 = parse_param(options,'beta2',0.1);
% tau = beta2/beta;

beta = beta * (b/B);
% options.tau = tau;

gamma = parse_param(options,'gamma',0.5);%weight
gamma = max(gamma,1e-4);
gamma = 1/(N*B*(1-gamma)/(N1*b*gamma) + 1);

%% Do fusion
disp('Start image fusion');

F = S*H;
X1 = X - ones(N1,1) * h0';
Z = gamma*G'*Y + (1-gamma)*X1*F';

options.regularization = regularization;
L = calc_regularization(reshape(X1,[rows1,cols1,b]), options);

options = insert_into_options(options,'g',g, ...
    'beta',beta,'gamma',gamma,'G',G,'Y',Y,'F',F,'X1',X1,'Z',Z, ...
    'extra',extra,'s',s,'L',L);

% R = optimize_gradient_descent(I, I1, options);
% R = fusion_with_couple_Local_LMM(I,I1,options);

if strcmp(regularization, 'LLMM')
    if 1 % reconstruct full reflectance data
        se_A = gamma*(G'*G) + beta*L;
        se_B = (1-gamma)*(F*F');
        se_C = Z;
        R = sylvester_krylov(se_A, se_B, se_C, sylvester_solver_mode);
    else % reconstruct abundances
        num_endm = vd(Y',5*10^-2); % M can be automatically defined, for example, by VD
        [~, indicies] = vca(Y', num_endm);
        M = Y(indicies,:);
        [U,Lambda] = eig(M*M');
        tmp = M*(1-gamma)*F*F'*M';
        [V,S] = eig((tmp + tmp')/2);
        Lambda1 = diag(sqrt(1./diag(Lambda)));
        se_A = gamma*(G'*G) + beta*L;
        se_B = Lambda1*U'*V*S*V'*U*Lambda1;
        se_B = (se_B + se_B) / 2;
        se_C = Z*M'*U*Lambda1;
        Z = sylvester_krylov(se_A, se_B, se_C, 'adaptive');
        A = Z*Lambda1*U';
        R = A*M;
    end
elseif strcmp(regularization, 'LMM')
    kappa = parse_param(options,'kappa',0.1);
    kappa = kappa * (b/B);
    options.kappa = kappa;
    
    R = fusion_with_LMM_constraint(I, I1, options);
end

%% pad bad bands
if isfield(options,'bbl') && ~all(options.bbl) 
    R1 = nan(size(R,1),length(options.bbl));
    R1(:,options.bbl) = R;
    R = R1;
end

I2 = reshape(R, [rows1,cols1,size(R,2)]);

if pad_color
    I2 = remove_rgb_margin(I2, extra);
end

t_elapsed = toc(t_start);
disp(['Elapsed time for image fusion is ', num2str(t_elapsed)]);
end



