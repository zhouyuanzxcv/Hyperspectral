function [X,A,row,col] = reshape_hsi(I,A)
%RESHAPE_HSI Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    A = [];
end

if ndims(I) == 3
    [row,col,B] = size(I);
    N = row*col;
    if ~isempty(A)
        [row,col,M] = size(A);
        A = reshape(A, N, M);
    end
    X = reshape(I, N, B);
else
    [N,B] = size(I);
    row = N;
    col = 1;
    X = I;
end

