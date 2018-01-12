function [X_p] = gmm_project(X, mapping)
[N,B] = size(X);
X_p = (X - repmat(mapping.mean,N,1)) * mapping.M;
end