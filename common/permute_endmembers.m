function [P,M2_1] = permute_endmembers(M1,M2)
% PERMUTE_ENDMEMBERS
% Permute the rows of M2 such that the difference between M1 and M2 is
% minimized.
% Input:
%   M1 - ground truth endmembers
%   M2 - algorithm endmembers
% Output:
%   P - size(M1,1) by size(M2,1) permutation matrix
%   M2_1 - permuted M2
P = zeros(size(M1,1),size(M2,1));

if size(M1,1) < 10 && size(M2,1) < 10 && size(M1,1) == size(M2,1)
    % use brutal force
    best_err = Inf;
    best_p = [];

    P1 = perms([1:size(M1,1)]);
    if size(M1,1) > size(M2,1)
        P1 = unique(P1(:,1:size(M2,1)),'rows');
    end

    for i = 1:size(P1,1)
        p = P1(i,:);
        M1_1 = M1(p,:);
        err = calc_error(M1_1,M2);
        if err < best_err
            best_err = err;
            best_p = p;
        end
    end
    for i = 1:length(best_p)
        P(best_p(i),i) = 1;
    end
else
    % Use the greedy algorithm
    M1_1 = M1;
    M2_1 = M2;
    while ~isempty(M2_1) && ~isempty(M1_1)
        errors = [];
        indices = [];
        for i = 1:size(M2_1,1)
            m2 = M2_1(i,:);
            [C,ind] = min(nanmean(abs(repmat(m2,size(M1_1,1),1) - M1_1),2));
            errors(i) = C;
            indices(i) = ind;
        end
        [~,ind2] = min(errors);
        ind1 = indices(ind2);
        % find indices in the original matrix
        [~,ind1_ori] = min(nanmean(abs(repmat(M1_1(ind1,:),size(M1,1),1) - M1),2));
        [~,ind2_ori] = min(nanmean(abs(repmat(M2_1(ind2,:),size(M2,1),1) - M2),2));
        
        P(ind1_ori,ind2_ori) = 1;
        M1_1(ind1,:) = [];
        M2_1(ind2,:) = [];
    end
end

for i = 1:size(P,1)
    if all(P(i,:)==0)
        P(i,:) = NaN;
    end
end
        

M2_1 = P*M2;
for j = 1:size(M2_1,1)
    if all(isnan(M2_1(j,:)))
        M2_1(j,:) = -1;
    end
end


function err = calc_error(M1,M2)
err = nanmean(nanmean(abs(M1-M2)));