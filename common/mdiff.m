function [abs_err, rel_err] = mdiff(A, B, verbose)
%MDIFF Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    verbose = 1;
end

if iscell(A)
    e1s = zeros(1,length(A));
    e2s = e1s;
    for i = 1:length(A)
        [e1,e2] = mdiff_i(A{i}, B{i});
        e1s(i) = e1;
        e2s(i) = e2;
    end
    e1 = mean(e1s);
    e2 = mean(e2s);
else
    [e1,e2] = mdiff_i(A,B);
end

abs_err = e1;
rel_err = e1/e2;

if verbose
    disp(['Absolute: ', num2str(abs_err), ', Relative: ', num2str(rel_err)]);
end

end

function [e1,e2] = mdiff_i(A, B)
e = abs(A - B);
e1 = mean(e(:));

e2 = mean(abs(A(:)));
end