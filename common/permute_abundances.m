function P = permute_abundances(A1,A2)
%PERMUTE_ABUNDANCES Permute abundance maps in A2 such that they match the
%abundance maps in A1. Warning: it assumes that most elements in the
%abundance maps are 1s. To permute arbitrary values, use
%permute_endmembers.
% Input:
%   A1: N by M_1 matrix of the abundance maps to be matched.
%   A2: N by M_2 matrix of the abundance maps to be permuted.
% Note: M_1 = M_2 or M_1 < M_2 or M_1 > M_2
% Output:
%   P: Permutation matrix
A1_1 = A1;
A2_1 = A2;
P = zeros(size(A1,2),size(A2,2));
while ~isempty(A1_1) && ~isempty(A2_1) && ~all(all(isnan(A2_1)))
    closeness = [];
    indices = [];
    for i = 1:size(A2_1,2)
        A = A2_1(:,i);
        [ind,close] = find_closest_abundance_map(A1_1,A);
        closeness(i) = close;
        indices(i) = ind;
    end
    [~,ind_2] = max(closeness);
    ind_1 = indices(ind_2);
    A_1 = A1_1(:,ind_1);
    A_2 = A2_1(:,ind_2);
    
    [~,ind1] = min(mean(abs(repmat(A_1,1,size(A1,2)) - A1),1));
    [~,ind2] = min(mean(abs(repmat(A_2,1,size(A2,2)) - A2),1));
    P(ind1,ind2) = 1;
    
    A1_1(:,ind_1) = [];
    A2_1(:,ind_2) = [];
end

for i = 1:size(P,1)
    if all(P(i,:)==0)
        P(i,:) = NaN;
    end
end


function [ind,max_closeness] = find_closest_abundance_map(A1,A)
for i = 1:size(A1,2)
    closeness(i) = mean(A(A1(:,i)==1));
end
[max_closeness,ind] = max(closeness);
