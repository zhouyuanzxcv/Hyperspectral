function [ output_args ] = test_AAM( input_args )
%TEST_AAM Summary of this function goes here
%   Detailed explanation goes here
%TEST_MESMA Summary of this function goes here
%   Detailed explanation goes here
addpath('../GMM');

close all;

dataset = '7';
[endmembers,I,Y,R_gt,A_gt,names,wl] = prepare_supervised_unmixing(dataset);
[Y,~,rows,cols] = reshape_hsi(I);

M = length(endmembers);

for j = 1:M
    endmembers{j} = endmembers{j}';
end

t_start = tic;
if 0 % Use AAM
    [index1, A, reconstruction1, error1] = AAM(Y', endmembers);
else % Use MESMA
    [index2, A, reconstruction2, error2] = MESMA_bruteforce(Y', endmembers);
end
t_elapsed = toc(t_start);
disp(['Elapsed time for AAM or MESMA is ', num2str(t_elapsed), ' seconds']);

save('result_AAM.mat','A');

mdiff(A,A_gt);
show_abundances(A,rows,cols);

end


