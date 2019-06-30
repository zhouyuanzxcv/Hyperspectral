function [U,V] = optimize_wrt_V(rgb1, I1, wl, options)
%OPTIMIZE_WRT_V Summary of this function goes here
%   Detailed explanation goes here
degree = options.degree;
t = options.t;
s = options.s;
T = create_T(degree, t);
sigma = options.sigma;

I2 = double(rgb1) / 255;
I2 = transform(I2,T,s,sigma,size(I1,2),size(I1,1),options);

I_sel = select_relevant_bands(I1,wl,options.srf_range);

% use replicate for the boundary condition of I
Ix = zeros(size(I_sel));
Iy = zeros(size(I_sel));
for i = 1:size(I_sel,3)
    [fx,fy] = calc_gradient(I_sel(:,:,i));
    Ix(:,:,i) = fx;
    Iy(:,:,i) = fy;
end

U = options.U;
V = options.V;

xmin = 0.5;
xmax = size(I_sel,2) - 0.5;
ymin = 0.5;
ymax = size(I_sel,1) - 0.5;
x = linspace(xmin, xmax, size(I_sel,2));
y = linspace(ymin, ymax, size(I_sel,1));
[X,Y] = meshgrid(x,y);

I_i = zeros(size(I_sel));
Ix_i = zeros(size(I_sel));
Iy_i = zeros(size(I_sel));
I_i_h = I_i;
Ix_i_h = Ix_i;
Iy_i_h = Iy_i;

der_Us = zeros(size(I2));
der_Vs = zeros(size(I2));

old_U = U;
old_V = V;

delta_t = parse_param(options, 'optimize_V_delta_t', 1);
% alpha controls the smoothness of v
alpha = parse_param(options, 'optimize_V_alpha', 0.05);

% rescale alpha since it is set based on 3 bands
alpha = (alpha / 3) * size(I2,3); 

beta = 0; % beta controls the magnitude of v
iter_num = 500;
errors = nan(1,iter_num);
for iter = 1:iter_num
    XI = X + U;
    YI = Y + V;
    
    % interpolate the hyperspectral image
    I_i = interp2_hyper_internal(I_sel,X,Y,XI,YI,[xmin,xmax],[ymin,ymax]);
    
    % update the SRF
    if mod(iter,100) == 1 % make sure H1 is constructed for the first iter
        I3 = cat(3, ones(size(I_i,1), size(I_i,2)), I_i);

        H1 = solve_for_H1(reshape_hsi(I2), reshape_hsi(I3), options);
        H = H1(2:end,:);
    end

    errors(iter) = calc_err(I_i, I2, H1);

    % interpolate the gradient
    for i = 1:size(I_sel, 3)
%         I_i(:,:,i) = interp2(X, Y, I_sel(:,:,i), XI, YI, 'spline');
        Ix_i(:,:,i) = interp2(X, Y, Ix(:,:,i), XI, YI, 'linear', 0);
        Iy_i(:,:,i) = interp2(X, Y, Iy(:,:,i), XI, YI, 'linear', 0);
    end
    
    for j = 1:size(I2,3)
        for i = 1:size(I_sel,3)
            I_i_h(:,:,i) = H(i,j) * I_i(:,:,i);
            Ix_i_h(:,:,i) = H(i,j) * Ix_i(:,:,i);
            Iy_i_h(:,:,i) = H(i,j) * Iy_i(:,:,i);
        end
        left_multiple = H1(1,j) + sum(I_i_h, 3) - I2(:,:,j);
        der_Us(:,:,j) = left_multiple .* sum(Ix_i_h, 3);
        der_Vs(:,:,j) = left_multiple .* sum(Iy_i_h, 3);
    end
    
    % use symmetric for the boundary condition of v
    der_U = 2 * sum(der_Us, 3) - 2*alpha * calc_del2(U) + 2*beta * U;
    der_V = 2 * sum(der_Vs, 3) - 2*alpha * calc_del2(V) + 2*beta * V;
    U = U - delta_t * der_U;
    V = V - delta_t * der_V;
    
    
    % A strict threshold 1e-3 is required to get good results for the
    % simulated Pavia dataset
    if max(abs(U(:) - old_U(:))) < 1e-4 && max(abs(V(:) - old_V(:))) < 1e-4
        disp(['Estimating V stopped at iteration ',num2str(iter)]);
        break
    end
        
    old_U = U;
    old_V = V;
end

if iter >= iter_num
    disp('Estimating V runs out of iterations');
end

function [Ix,Iy] = calc_gradient(I)
[Ix,Iy] = gradient(I);

% use replicate for the boundary condition of I
Ix(:,1) = Ix(:,1)/2;
Ix(:,end) = Ix(:,end)/2;

Iy(1,:) = Iy(1,:)/2;
Iy(end,:) = Iy(end,:)/2;

function delta_U = calc_del2(U)
U1 = append_boundary(U);
delta_U = del2(U1);
delta_U = delta_U(2:end-1,2:end-1);

function f1 = append_boundary(f)
[rows,cols] = size(f);
f1 = zeros(rows+2,cols+2);
f1(2:end-1,2:end-1) = f;
f1(2:end-1,1) = f(:,2);
f1(2:end-1,end)  = f(:,end-1);
f1(1,2:end-1) = f(2,:);
f1(end,2:end-1) = f(end-1,:);

function err = calc_err(I1, I2, H1)
N = size(I1,1) * size(I1,2);
Y1 = [ones(N,1), reshape_hsi(I1)];
X = reshape_hsi(I2);
err = sum(sum((Y1 * H1 - X).^2));

function I1 = interp2_hyper_internal(I,X,Y,XI,YI,xrange,yrange)
%INTERP2_HYPER Summary of this function goes here
%   Detailed explanation goes here
% xrange = [min(X(:)), max(X(:))]
% yrange = [min(Y(:)), max(Y(:))]

I1 = zeros(size(XI,1),size(XI,2),size(I,3));

% use replicate for the boundary condition of I
XI(XI < xrange(1)) = xrange(1);
XI(XI > xrange(2)) = xrange(2);
YI(YI < yrange(1)) = yrange(1);
YI(YI > yrange(2)) = yrange(2);

% this part is used for nonrigid transformation of hyperspectral images.
for i = 1:size(I,3)
    I1(:,:,i) = interp2(X, Y, I(:,:,i), XI, YI, 'linear');
end

