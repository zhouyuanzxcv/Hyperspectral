function R = restore_endm_by_bbl(R, bbl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if iscell(R)
    for i = 1:length(R)
        R{i} = restore_from_bbl(R{i}, bbl);
    end
else
    R = restore_from_bbl(R, bbl);
end

end

function R1 = restore_from_bbl(R, bbl)
[M,B1] = size(R);
R1 = NaN(M, length(bbl));
R1(:,bbl) = R;

end