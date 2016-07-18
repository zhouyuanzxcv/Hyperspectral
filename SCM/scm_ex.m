function [sigma,mu,var_dirs,var_amts] = scm_ex(Y,A,R,sigma0,mu0,sigma_max)
%SCM_EX Summary of this function goes here
%   Detailed explanation goes here
[N,M] = size(A);
[~,B] = size(R);
sigma_min = 1e-6;

C = kron((A'*A),eye(B));
S = (mu0^2/(sigma0^2))*eye(M*B,M*B);
    
h = (Y-A*R)'*A;
h = h(:);
gamma = 1/(mu0^2);
lsq = sum(sum((Y-A*R).^2));
gamma0Inv = 1/(N*B)*lsq;
obj_fun_sigma = @(S,gamma) gamma*lsq - gamma*h'*((S+C)\h) ...
    + logdet(S+C) - logdet(S) - N*B*log(gamma);
calc_ignores = @(S,gamma) [(gamma*h'*((S+C)\h)) / (gamma*lsq), ...
    (logdet(S+C) - logdet(S)) / (gamma*lsq)];
disp(['Inital ignores are: ',num2str(calc_ignores(S,gamma))]);
if 0
    % U1'*N1 has a much smaller frobenius norm than that of N1
    global A_gt; global R_gt;
    N = Y-A*R; N1 = Y-A_gt*R_gt;
    [U,D,V] = svd(A,0); [U1,D1,V1] = svd(A_gt,0);
    sum(sum(N.^2)), sum(sum((U'*N).^2))
    sum(sum(N1.^2)), sum(sum((U1'*N1).^2))
end
Es = obj_fun_sigma(S,gamma);
delta_t0 = 1e-9;

T = S + C;
Th = T\h;

disp('Start estimating the uncertainty...');
for i = 1:200
    % update S
    S2_old = S;

    GammaTinv = gamma*Th*Th';
    Tinv = inv(T);
    E0 = obj_fun_sigma(S,gamma);
    delta_t = delta_t0;

    der = zeros(M*B,M*B);
    for j = 1:M     
        inds = (j-1)*B+1:j*B;
        Sj = S(inds,inds);
        der(inds,inds) = GammaTinv(inds,inds) - inv(Sj) + Tinv(inds,inds);
    end


    while 1
        S2_new = S - delta_t*der;
        S2_new = (S2_new + S2_new')/2;
        for j = 1:M
            inds = (j-1)*B+1:j*B;
            S3 = S2_new(inds,inds);

            [V,D] = eig(S3);
            d = diag(D);

            epsilon = 1/(gamma*sigma_max^2);
            infinity = 1/(gamma*sigma_min^2);
            d(d<epsilon) = epsilon;
%                 d(d>infinity) = infinity;

            D = diag(d);
            S2_new(inds,inds) = V*D*V';
        end

        E0 = [E0; obj_fun_sigma(S2_new,gamma)];
        if E0(end) >= E0(end-1)
            S = S2_old;
            break;
        else
            S2_old = S2_new;
            delta_t = delta_t*10;
        end
    end
    % Update gamma
    T = S + C;
    Th = T\h;
    gamma = gamma0Inv - 1/(N*B)*h'*Th;
    gamma = 1/gamma;
    % Check the amounts of ignored terms for estimating A and R
    if 1
        total_amt = gamma*lsq;
        amt1 = gamma*h'*((S+C)\h);
        amt2 = logdet(S+C) - logdet(S);
        ratio1(i) = amt1/total_amt;
        ratio2(i) = amt2/total_amt;
    end
    % Check convergence
    Es = [Es;obj_fun_sigma(S,gamma)];  
    if Es(end-1) - Es(end) <= 1e-6*(Es(1) - Es(2))
        disp(['Finishes at iteration ',num2str(i),', current energy: ',num2str(Es(end))]);
        disp(['Final ignores are: ',num2str(calc_ignores(S,gamma))]);
        break;
    end

    if mod(i,5) == 0
        disp(['Process iteration ',num2str(i),', current energy: ',num2str(Es(end))]);
    end
end

mu = gamma^(-1/2);
d = (1/mu) * ones(B,1);

[sigma,var_dirs,var_amts] = calc_uncertainty_range(S,d);

