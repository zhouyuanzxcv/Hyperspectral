function A2 = project_to_simplex(A)
%PROJECT_TO_SIMPLEX Summary of this function goes here
%   Detailed explanation goes here
[N,M] = size(A);
[A1,IX] = sort(A',1,'descend');
J = ones(N,1)*(1:M); J = J';
A1cum = cumsum(A1,1)-1*ones(M,N);
A1minus = (A1 - (A1cum ./ J)) > 0;

rho = max(double(A1minus).*J,[],1);
theta = A1cum(sub2ind([M N],rho,(1:N)))./rho;

A1 = A1 - ones(M,1)*theta;
A1(A1<0) = 0;

% reorder
A2 = zeros(N,M);
for i = 1:M
    A2(:,i) = A1(IX==i);
end

