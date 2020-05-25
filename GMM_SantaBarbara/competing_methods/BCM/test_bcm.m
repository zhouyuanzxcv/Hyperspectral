function [ output_args ] = test_bcm(dataset)
%TEST_BCM Summary of this function goes here
%   Detailed explanation goes here
addpath('../GMM');
if nargin < 1
    dataset = '001';
end
    
options = struct('force_positive_I',1);
[endmembers,I,Y,R_gt,A_gt,names,wl] = prepare_supervised_unmixing(dataset,options);

M = size(A_gt,2);
[rows,cols,B] = size(I);

for j = 1:M
    endmembers{j} = endmembers{j}';
end

[Parameters] = BCMParameters(endmembers);

disp('BCM-Spectral-QP Unmixing...');

t_start = tic;
[A] = BCM(I, Parameters, 1);%BCM-Spectral-QP
t_elapsed = toc(t_start);
disp(['BCM lasts ', num2str(t_elapsed), ' seconds']);

R = calc_endmembers(endmembers);

save(['result_bcm_spectral_qp_',dataset,'.mat'],'A','R');

mdiff(A,A_gt);

show_abundances(A,rows,cols,'BCM-Spectral-QP');

if 0
    disp('BCM-Spatial-QP Unmixing...');
    [A] = BCM(I, Parameters, 3);%BCM-Spatial-QP
    R = calc_endmembers(Y,A);
    show_abundances(A,rows,cols,'BCM-Spatial-QP');

    save('result_bcm_spatial_qp.mat','A','R');
end

function R = calc_endmembers(endmembers)
M = length(endmembers);
B = size(endmembers{1},1);
R = zeros(M,B);
for j = 1:M
    R(j,:) = mean(endmembers{j},2)';
end


