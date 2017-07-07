function [error_avg,error_A] = calc_abundance_error(A1,A2)
%CALC_ABUNDANCE_ERROR Summary of this function goes here
%   Detailed explanation goes here
if ~isfloat(A1)
    A1 = double(A1);
end

if length(unique(A1(:))) == 2
    is_real_dataset = 1;
else
    is_real_dataset = 0;
end

P = permute_abundances(A1,A2);
if nanmean(nanmean(abs(eye(size(A1,2)) - P))) ~= 0
    A2_1 = (P*A2')';
else
    A2_1 = A2;
end

if is_real_dataset
    inds = any(A1==1,2);
    error_A = nanmean(abs(A1(inds,:) - A2_1(inds,:)));
else
    error_A = nanmean(abs(A1-A2_1));
end

error_avg = nanmean(error_A);