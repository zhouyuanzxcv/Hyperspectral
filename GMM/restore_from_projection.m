function [mus,sigmas,R] = restore_from_projection(mus_p,sigmas_p,R_p,c,E)
% restore from projection.
% input:
%   mus_p - 1 by M cell array (each cell is a K_j by B matrix)
%   sigmas_p - 1 by M cell array (each cell is a B by B by K_j matrix)
%   R_p - M by d projected data or M by d by N projected data
%   c - 1 by B center of the PCA
%   E - B by d projection matrix
% output:
%   mus, sigmas, R - restored versions of mus_p, sigmas_p and R_p

if ~isempty(mus_p) && ~isempty(sigmas_p)
    M = length(mus_p);
elseif ~isempty(R_p)
    M = size(R_p,1);
end
mus = cell(1,M);
sigmas = cell(1,M);
B = size(E,1);
R = zeros(M,B);

if ~isempty(mus_p) && ~isempty(sigmas_p)
    for j = 1:M
        K_j = size(mus_p{j},1);
        mus{j} = zeros(K_j,B);
        sigmas{j} = zeros(B,B,K_j);

        for k = 1:size(mus_p{j},1)
    %         mus{j}(k,:) = solve_linsys_mu(E,mus_p{j}(k,:)')' + c;
    %         sigmas{j}(:,:,k) = solve_linsys_sigma(E, sigmas_p{j}(:,:,k));
            mus{j}(k,:) = recover_mu(E,c,mus_p{j}(k,:));
            sigmas{j}(:,:,k) = recover_sigma(E,sigmas_p{j}(:,:,k));
        end
    end
end

if ~isempty(R_p)
    if ndims(R_p) == 2
        R = recover_mu(E,c,R_p);
    elseif ndims(R_p) == 3
        R = zeros(M,B,size(R_p,3));
        for i = 1:size(R_p,3)
            R(:,:,i) = recover_mu(E,c,R_p(:,:,i));
        end
    end
end

function mu = recover_mu(E,c,mu0)
mu = E*mu0' + repmat(c',[1, size(mu0,1)]);
mu = mu';

function sigma = recover_sigma(E, sigma0)
sigma = E*sigma0*E' + 1e-9*eye(size(E,1));

function x = solve_linsys_mu(E,mu0)
EE1 = E*E';
B = size(EE1,1);
EE1 = EE1 + 1e-9*eye(B);
x = EE1 \ (E* mu0);

function sigma = solve_linsys_sigma(E,sigma0)
[B,d] = size(E);
w = 1e-9;
EE = kron(E,E);
tmp1 = EE * sigma0(:);
tmp2 = (eye(d^2) + (1/w)*EE'*EE) \ (EE' * tmp1);
sigma = (1/w) * tmp1 - (1/w^2) * EE * tmp2;
sigma = reshape(sigma,B,B);

