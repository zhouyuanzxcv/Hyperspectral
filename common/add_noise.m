function I1 = add_noise(I, Y_noise)
% Input:
%   I - rows by cols by B image data
%   Y_noise - scalar or B by 1 vector of standard deviations
% Output:
%   I1 - noised contaminated image
[rows,cols,B] = size(I);
N = rows * cols;
if length(Y_noise) == 1
    Y_noise = repmat(Y_noise, [B,1]);
end

X = reshape(I, [rows*cols, B]);
Y = X + randn(N,B) .* repmat(Y_noise(:)', [N 1]);
I1 = reshape(Y, [rows,cols,B]);
end

