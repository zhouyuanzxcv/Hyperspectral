function D = calc_distance_matrix(A,B)
%CALC_DISTANCE Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   A: N by D data
%   B: M by D data
% Output:
%   D: N by M distance matrix
if size(A,2) ~= size(B,2)
    disp('Error: the dimensions do not match.');
end

A1 = sum(A.^2,2); 
B1 = sum(B.^2,2); 
AB = A*B'; 
D = sqrt(repmat(A1,[1 length(B1)]) + repmat(B1',[length(A1) 1]) - 2*AB);

