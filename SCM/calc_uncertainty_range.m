function [sigma,var_dirs,var_amts] = calc_uncertainty_range(S,d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
MB = size(S,1);
B = length(d);
M = MB/B;

sigma = zeros(B,B,M); % covariance matrix
var_dirs = zeros(M,B);
var_amts = zeros(M,1);

for j = 1:M
    inds = (j-1)*B+1:j*B;
    sigma1 = inv(diag(d)*S(inds,inds)*diag(d));
    sigma(:,:,j) = (sigma1 + sigma1')/2; %It is symmetric theoretically

    % computer the major variations
    sigma1 = sigma(:,:,j);
    [V,D] = eig(sigma1);
    d1 = diag(D);
    [c,ind] = max(d1);
    var_dir = V(:,ind);
    var_amt = sqrt(d1(ind));
    var_dirs(j,:) = var_dir';
    var_amts(j) = var_amt;
end

