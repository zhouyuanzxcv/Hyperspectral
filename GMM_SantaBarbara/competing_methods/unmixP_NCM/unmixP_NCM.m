function [P] = unmixP_NCM(X,E,Sigma,Parameters)
%% Input:
%        X - d-by-N HSI data, matrix
%        E - d-by-M endmember mean set, matrix
%        Sigma - d-by-d-by-M endmember covariance set, matrix
%  Output:
%        P - M-by-N proportion matrix
% Author: Alina Zare et al. Rewriten by Sheng Zou
% Department of Electrical and Computer Engineering, University of Florida
% 09/07/2017

%% Initialization of P
[d,N] = size(X);
M = size(E,2);
P = 1/M.*ones(M,N); % intialize proportion values to be 1/M

%Initialize all Likelihood and Prior Values
LogLikelihoodOld          = ComputeLogLikelihoodAll(X, E, P, Sigma, N);


for iteration = 2:Parameters.NumberIterations+1
    Y = randg(1, N, M) ;
    v = sum(Y,2);
    samples = (Y./repmat(v,[1,size(Y,2)]))';
    LogLikelihoodNew      = ComputeLogLikelihoodAll(X, E, samples, Sigma, N);
    Ratios                = exp(LogLikelihoodNew - LogLikelihoodOld);
    rands                 = rand(1,N);
    Vals                  = rands < Ratios;
    Vrep                  = repmat(Vals,M,1);
    P                     = samples.*Vrep  + P.*(1-Vrep);
    
    LogLikelihoodOld = (1-Vals).*LogLikelihoodOld + (Vals).*LogLikelihoodNew;
    if mod(iteration, round(Parameters.NumberIterations/10)) == 0
        disp(strcat('iteration = ',num2str(iteration)));
        disp(strcat('loglikelihood = ',num2str(sum(LogLikelihoodOld))));
    end
end

end

function [LogLikelihoodAll] = ComputeLogLikelihoodAll(X, E, P, cov, N)
% compute log likelihood of all points

term1 = zeros(1,N);
term2 = zeros(1,N);

for s = 1:size(E,2)
    statement(s)=isdiag(squeeze(cov(:,:,s))) && length(unique(diag(cov(:,:,s))))==1; % check if all the covariance matrices are diagonal and isotropic
end

if mean(statement) ==1 % all diagonal and isotropic
    for z = 1:size(P,1)
        a(z) = unique(diag(cov(:,:,z)));
    end
    term1 = -.5*size(X,1)*log(a*(P.^2));
    term2 = -.5*(sum((X-E*P).^2)./(a*(P.^2)));
else
    
    parfor t = 1:N
        P3D = zeros(1,1,size(E,2));
        P3D(:,:,1:size(E,2)) = P(:,t).^2;
        P3D_full = repmat(P3D,size(X,1),size(X,1));
        term1(t) = -.5*logdet(sum(P3D_full.*cov,3));
        term2(t) = -.5*(X(:,t)-E*P(:,t))'*(sum(P3D_full.*cov,3)\(X(:,t)-E*P(:,t)));
    end
    
end

LogLikelihoodAll = term1 + term2;
end



