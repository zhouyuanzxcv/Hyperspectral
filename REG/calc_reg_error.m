function [mean_val,translations] = calc_reg_error(rgb1,I1,T,s,U,V,T2,s2,U2,V2)
%CALC_REG_ERROR Summary of this function goes here
%   Detailed explanation goes here
[rows,cols,B] = size(I1);
[rows1,cols1,b] = size(rgb1);


% The two approaches give almost the same error
if 1
    % calculate centers in the hi-resolution domain based on ground truth
    centers = calc_centers(rows, cols, T, s, U, V);

    % transform back to the hyperspectral domain using estimated transformation
    centers1 = transform_back(centers, rows, cols, T2, s2, U2, V2);

    % calculate expected coordinates
    [X,Y] = meshgrid((0.5:cols-0.5),(0.5:rows-0.5));
    centers2 = cat(3,X,Y);

    % compare transformed to expected
    translations = sqrt(sum((centers1 - centers2).^2, 3));
else
    % calculate centers in the hi-resolution domain based on ground truth
    centers = calc_centers(rows, cols, T, s, U, V);
    
    % calculate centers in the hi-resolution domain based on estimated
    centers1 = calc_centers(rows, cols, T2, s2, U2, V2);
    
    % calculate distance and scale it to the hyperspectral domain
    translations = sqrt(sum((centers1 - centers).^2,3))/mean(s);
end

translations = translations(:);
mean_val = mean(translations);

end

function centers1 = transform_back(centers, rows, cols, T, s, U, V)
centers1 = inv(T) * cat(1, reshape(centers, [rows*cols,2])', ones(1,rows*cols));
centers1(1,:) = centers1(1,:) / s(1);
centers1(2,:) = centers1(2,:) / s(2);
centers1(3,:) = [];
centers1 = reshape(centers1', [rows,cols,2]);
centers1(:,:,1) = centers1(:,:,1) + U;
centers1(:,:,2) = centers1(:,:,2) + V;
end

function centers = calc_centers(rows, cols, T, s, U, V)
x = (0.5:cols-0.5);
y = (0.5:rows-0.5);
[X,Y] = meshgrid(x,y);
X = X - U;
Y = Y - V;
X = X * s(1);
Y = Y * s(2);
coord = cat(3,X,Y);
centers = T * cat(1, reshape(coord, [rows*cols, 2])', ones(1,rows*cols));
centers(3,:) = [];
centers = reshape(centers', [rows,cols,2]);
end

