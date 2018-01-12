function [A,J2,A2] = create_adjacency_graph(X,method,param,force_sym)
%CREATE_ADJACENCY_GRAPH Summary of this function goes here
%   Detailed explanation goes here
if strcmp(method,'epsball')
    fh = @(Z) find_indices_by_epsilon(Z,param);
    [J2,A2] = create_adjacency_graph_internal(X,fh);
elseif strcmp(method,'nn')
    if 0 % my implementation
        fh = @(Z) find_indices_by_nn(Z,param);
        [J2,A2] = create_adjacency_graph_internal(X,fh);
    else % Matlab implementation is faster
        [D,I] = pdist2(X,X,'euclidean','Smallest',param+1);
        D(1,:) = [];
        I(1,:) = [];
        J2 = cell(size(D,2),1);
        A2 = cell(size(D,2),1);
        for i = 1:length(J2)
            J2{i} = I(:,i)';
            A2{i} = D(:,i)';
        end
    end
else
    disp('The method parameter can only take epsball or nn');
    A = []; J2 = []; A2 = [];
    return;
end
A = cellarr2sparse(J2,A2,length(J2),length(J2));
if force_sym
    A = max(A,A');
end

function [J2,A2] = create_adjacency_graph_internal(X,fh_find_indices)
step = 100;
n = size(X,1);
J2 = cell(n,1);
A2 = cell(n,1);
for i1 = 1:step:n
    i2 = i1 + step - 1;
    if (i2 > n)
        i2 = n;
    end;

    XX = X(i1:i2,:);
    D = calc_distance_matrix(XX,X);
    [Z,I] = sort(D,2);

    Z = Z(:,2:end); % remove itself as closest neighbor
    I = I(:,2:end);
    for i = i1:i2
        ind = fh_find_indices(Z(i-i1+1,:));
%         ind = find(Z(i-i1+1,:) < epsilon);
        jj = I(i-i1+1,ind);
        Z1 = Z(i-i1+1,ind);
        J2{i} = jj;
        A2{i} = Z1;
    end
end

function ind = find_indices_by_epsilon(Z,epsilon)
ind = find(Z < epsilon);

function ind = find_indices_by_nn(Z,nn)
ind = [1:min(nn,length(Z))];
