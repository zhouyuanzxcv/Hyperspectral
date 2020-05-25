function [average_mad, average_cc, extra] = calc_error_measure(A_gt_total, As)
%CALC_ERROR_MEASURE Summary of this function goes here
%   Detailed explanation goes here
A_total = zeros(length(As), size(As{1},2));
for i = 1:length(As)
    A = As{i};
    A_total(i,:) = mean(A,1);    
end

disp('MAD:');
errors = abs(A_gt_total - A_total);
endm_error_mean = mean(errors,1);
endm_error_std = std(errors);
endm_error_rmse = sqrt(mean(errors.^2, 1));
disp(['Mean: ',num2str(endm_error_mean)]);
disp(['Std: ',num2str(endm_error_std)]);
average_mad = mean(errors(:));
disp(['Total average: ',num2str(average_mad)]);

disp('R.^2:');
R = zeros(1,size(A_gt_total,2));
for i = 1:size(A_gt_total,2)
    R1 = corrcoef(A_gt_total(:,i), A_total(:,i));
    R(i) = R1(1,2);
end
R2 = R.^2;
disp(['Endmember: ',num2str(R2)]);
average_cc = mean(R2);
disp(['Average: ',num2str(average_cc)]);

extra = [];
extra.A_gt_total = A_gt_total;
extra.A_total = A_total;
extra.endm_error_mean = endm_error_mean;
extra.endm_error_std = endm_error_std;
extra.endm_error_rmse = endm_error_rmse;
extra.R2 = R2;

end

