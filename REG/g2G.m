function G = g2G(C,g,N)
% N = length(C);
% R2 = length(g);
% I = zeros(N*R2,1);
% J = zeros(N*R2,1);
% V = zeros(N*R2,1);
% for i = 1:N
%     inds = (i-1)*R2+1:i*R2;
%     [Ii,Ji,Vi] = find(g'*C{i});
%     I(inds) = i*ones(R2,1);
%     J(inds) = Ji;
%     V(inds) = Vi;
% end
% G = sparse(I,J,V,N,N*R2);

[R,N1N] = size(C);
N1 = N1N/N;
vecG1 = C'*sparse(g);
G1 = reshape(vecG1,N1,N);
G = G1';