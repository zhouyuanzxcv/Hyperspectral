function [d,S] = solve_for_d_S(Y,A,R,options)
[N,B] = size(Y);
M = size(A,2);

sigma0 = 0.1; % initial sigma
sigma_max = 1;
sigma_min = 1e-9;
delta_t0 = 1e-9;

disp('Start estimating D and S');
C = kron((A'*A),eye(B));
YAR = (Y-A*R)';
Z = YAR*A;
E = diag(sum(YAR.^2,2));
F1 = YAR*A;
F = zeros(M*B,B);
for i = 1:M
    F((i-1)*B+1:i*B,:) = diag(F1(:,i));
end

d = (1/N)*sum(YAR.^2,2);
d = 1./sqrt(d);

S = eye(M*B,M*B);
if 0
    u = 1/sqrt(B) * ones(B,1);
    S1 = sigma0^2*u*u' + 1e-4*eye(B);
    S1 = inv(S1);
else
    S1 = sigma0^(-2)*eye(B);
end

for j = 1:M
    inds = (j-1)*B+1:j*B;
    S(inds,inds) = diag(1./d) * S1 * diag(1./d);
end

errors2 = eval_obj_fun_D_S(S,d,C,YAR,Z,M,B,N);
for iter = 1:300
    % solve for S
    z = diag(d)*Z;
    z = z(:);
    if 0
        obj_fun_S = @(S) eval_obj_fun_S(S,C,z,M,B);
        fun_der_S = @(S) calc_der_S(S,C,z,M,B);
        proj_fun_S = @(S) project_to_spd(S,d,sigma_min,sigma_max,M,B);

        S = projected_gradient_descent(S, obj_fun_S, fun_der_S, ...
            proj_fun_S, delta_t0);
    else
        if 1
            proj_fun_S = @(S) project_to_spd(S,d,sigma_min,sigma_max,1,B);
            for j = 1:M
                obj_fun_Sj = @(Sj) eval_obj_fun_Sj(Sj,S,C,z,M,B,j);
                fun_der_Sj = @(Sj) calc_der_Sj(Sj,S,C,z,B,j);
                inds = (j-1)*B+1:j*B;
                Sj = S(inds,inds);
                Sj = projected_gradient_descent(Sj, obj_fun_Sj, fun_der_Sj, ...
                    proj_fun_S, 1e-5);
                S(inds,inds) = Sj;
            end
        else
            der_S = calc_der_S(S,C,z,M,B);
            proj_fun_S = @(S) project_to_spd(S,d,sigma_min,sigma_max,1,B);
            for j = 1:M
                obj_fun_Sj = @(Sj) eval_obj_fun_Sj(Sj,S,C,z,M,B,j);
                fun_der_Sj = @(Sj) calc_der_Sj_approx(Sj,der_S,B,j);
                inds = (j-1)*B+1:j*B;
                Sj = S(inds,inds);
                Sj = projected_gradient_descent(Sj, obj_fun_Sj, fun_der_Sj, ...
                    proj_fun_S, 1e-3);
                S(inds,inds) = Sj;
            end
        end
    end
    
    % solve for D
    Q = S + C;
    G = (1/N)*(E-F'*(Q\F));
    fcn_d = @(d) G*d - 1./d;
    fsolve_opts = optimset('Display','off');
    d = fsolve(fcn_d, 1./sqrt((1/N)*diag(E)), fsolve_opts);
        
    % calc error
    err = eval_obj_fun_D_S(S,d,C,YAR,Z,M,B,N);
    errors2 = [errors2;err];
    if test_convergence(errors2, 1e-8), break; end
    
    if mod(iter,10) == 0
        disp(['Process iteration ',num2str(iter),'. The objective ',...
            'function has value ',num2str(err)]);
    end
end


function S = project_to_spd(S,d,sigma_min,sigma_max,M,B)
touch_bd = 0;
S = (S + S')/2;
for j = 1:M
    inds = (j-1)*B+1:j*B;
    S3 = S(inds,inds);
    SigmaInv = diag(d) * S3 * diag(d);
    [V,D] = eig(SigmaInv);
    d1 = diag(D);
    low_bd = 1/(sigma_max^2);
    up_bd = 1/(sigma_min^2);
    
    if ~isempty(find(d1<low_bd, 1))
        touch_bd = touch_bd + 1;
%         break;
    end
    
    d1(d1<low_bd) = low_bd;
    d1(d1>up_bd) = up_bd;
    S(inds,inds) = diag(1./d)*V*diag(d1)*V'*diag(1./d);
end

if touch_bd >= 1
%     disp(['Uncertainty amount upper bound touched ',num2str(touch_bd),' times']);
    
    % Once uncertainty amount upper bound touched. It means the change of
    % the covariance is too large. Setting S to 0 will make the objective
    % function be Inf, thus revoke the delta_t.
    
%     if touch_bd >= 1
%         S = 0 * eye(M*B);
%     end
end


function val = eval_obj_fun_Sj(Sj,S,C,z,M,B,j)
inds = (j-1)*B+1:j*B;
S(inds,inds) = Sj;
val = eval_obj_fun_S(S,C,z,M,B);


function val = eval_obj_fun_S(S,C,z,M,B)
Q = S + C;
% R = chol(Q);
val1 = 0;
for j = 1:M     
    inds = (j-1)*B+1:j*B;
    val1 = val1 + logdet(S(inds,inds));
end

% val = -z'*(Q\z) + logdet(Q) - val1;
R = chol(Q);
y = R'\z;
val = -sum(y.^2) + 2*sum(log(diag(R))) - val1;

function val = eval_obj_fun_D_S(S,d,C,YAR,Z,M,B,N)
z = diag(d)*Z;
z = z(:);

val = eval_obj_fun_S(S,C,z,M,B);
val = val + sum(sum((diag(d)*YAR).^2)) - 2*N*sum(log(d));


function der = calc_der_S(S,C,z,M,B)
Q = S + C;
invQ = inv(Q);
y = invQ*z;
T = y*y';
der = zeros(M*B,M*B);
for j = 1:M
    inds = (j-1)*B+1:j*B;
    Sj = S(inds,inds);
    der(inds,inds) = T(inds,inds) + invQ(inds,inds) - inv(Sj);
end

% der_norm = norm(der,'fro');
% der = der/der_norm;


function der = calc_der_Sj(Sj,S,C,z,B,j)
inds = (j-1)*B+1:j*B;
S(inds,inds) = Sj;

Q = S + C;
invQ = inv(Q);
y = invQ*z;
T = y*y';
der = zeros(B,B);

der = T(inds,inds) + invQ(inds,inds) - inv(Sj);

der_norm = norm(der,'fro');
der = der/der_norm;

function der = calc_der_Sj_approx(Sj,der_S,B,j)
inds = (j-1)*B+1:j*B;
der = der_S(inds,inds);

der_norm = norm(der,'fro');
der = der/der_norm;

function d = solve_for_d1(Y,A,R,S)
[N,B] = size(Y);
M = size(A,2);

YAR = Y - A*R;
E = diag(sum(YAR.^2,1));

Q = S + kron((A'*A),eye(B));

F1 = A'*YAR;
F = zeros(M*B,B);
for i = 1:M
    F((i-1)*B+1:i*B,:) = diag(F1(i,:));
end
G = (1/N)*(E-F'*(Q\F));

fcn_d = @(d) G*d - 1./d;
d = fsolve(fcn_d, 1./sqrt((1/N)*diag(E)));