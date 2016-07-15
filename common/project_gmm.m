function [proj_mus, proj_sigmas] = project_gmm(mus, sigmas, mapping)
%PROJECT_GMM Summary of this function goes here
%   Detailed explanation goes here
M = length(mus);

proj_mus = cell(1,M);
proj_sigmas = cell(1,M);
for j = 1:M
    proj_mus{j} = (mus{j} - repmat(mapping.mean, size(mus{j},1), 1))*mapping.M;
    for k = 1:size(sigmas{j},3)
        sigma = mapping.M' * sigmas{j}(:,:,k) * mapping.M;
        proj_sigmas{j}(:,:,k) = (sigma + sigma') / 2; % force symmetry
    end
end


