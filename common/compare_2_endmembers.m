function [error_M,error_A,P] = compare_2_endmembers(M1,M2,A1,A2,m,n,...
    names,wl,show_figure,options)
%COMPARE_ENDMEMBERS Summary of this function goes here
%   Detailed explanation goes here
if nargin < 9
    show_figure = 0;
end

if nargin < 10
    options = [];
end

A1 = double(A1);

% use_uncertainty = 0;
var_dirs = [];
var_amts = [];

permute_criterion = 'auto';

if isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end

use_uncertainty = ~isempty(var_amts);

if isempty(M1) || isempty(A1) || size(A1,2) < size(A2,2)
    show_endmembers(M2,wl,names);
    show_abundances(A2,m,n);
    if use_uncertainty
        show_uncertainty_range([], M2, var_amts, var_dirs, wl, names, options);
    end
    
    error_M = [];
    error_A = [];
    P = eye(size(M2,1));
    return;
end


[M,B] = size(M1);

is_real_dataset = 0;
if length(unique(A1(:))) == 2
    is_real_dataset = 1;
end

if strcmp(permute_criterion, 'auto')
    if is_real_dataset
        permute_criterion = 'abundance';
    else
        permute_criterion = 'endmember';
    end
end

%% Permute the endmembers and abundances from the algorithm to accord with
% the ground truth ones.
if strcmp(permute_criterion, 'abundance')
    P = permute_abundances(A1,A2);
    M2_1 = P*M2;
elseif strcmp(permute_criterion, 'endmember')
    [P,M2_1] = permute_endmembers(M1,M2);
end
A2_1 = (P*A2')';


if use_uncertainty
    var_amts = P*var_amts;
    var_dirs = P*var_dirs;
end

%% Plot 2 groups of endmembers and show uncertainty range if possible
if show_figure
    show_uncertainty_range(M1, M2_1, var_amts, var_dirs, wl, names, options);
end

%% Calculate the error
error_A = calc_abundance_error(A1,A2_1);
error_M = nanmean(nanmean(abs(M1-M2_1)));

if show_figure
    disp(['Error for A: ',num2str(error_A)]);
    disp(['Error for M: ',num2str(error_M)]);
end

if show_figure
    show_abundances(A2_1,m,n);
end

%% If size(A1,2) < size(A2,2), make the permutation matrix P square
if size(P,1) < size(P,2)    
    missings = [];
    for i = 1:size(P,2)
        if all(P(:,i) == 0)
            missings = [missings,i];
        end
    end
    new_P1 = zeros(length(missings), size(P,2));
    for i = 1:length(missings)
        new_P1(i, missings(i)) = 1;
    end
    P = [P;new_P1];
end
