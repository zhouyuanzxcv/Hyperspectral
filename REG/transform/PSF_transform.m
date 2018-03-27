function I1 = PSF_transform(I, g, cols, rows, extra)
%PSF_TRANSFORM Summary of this function goes here
%   Detailed explanation goes here

% The input I is a padded image for the PSF
extra_x = extra(1);
extra_y = extra(2);

[rows1,cols1,b] = size(I);
[sy1,sx1] = size(g);
if sy1 == 1 && sx1 == 1 % the PSF is 1
    I1 = I;
    return;
end
    
sy = sy1-2*extra_y;
sx = sx1-2*extra_x;
assert(rows == (rows1-2*extra_y)/sy);
assert(cols == (cols1-2*extra_x)/sx);

N = rows*cols;
R2 = sy1*sx1;
% g1 = repmat(g(:), [1,b]);
g = g(:)';

if 1
    I1 = zeros(rows,cols,b);

%     for i = 1:N
%         [r,c] = ind2sub([rows,cols],i);
%         yi = (r-1)*sy + 1 : r*sy + 2*extra_y;
%         xi = (c-1)*sx + 1 : c*sx + 2*extra_x;
%         rect = I(yi,xi,:);
%         I1(r,c,:) = sum(g1 .* reshape(rect, [R2,b]), 1);
%     end
    
    % The code below is a translation of the above commented code but 10% faster.
    [r,c] = ind2sub([rows,cols],(1:N)');
    YI = repmat((r-1)*sy + 1, [1 sy + 2*extra_y]);
    YI = YI + repmat((0:sy + 2*extra_y - 1), [N 1]);
    XI = repmat((c-1)*sx + 1, [1 sx + 2*extra_x]);
    XI = XI + repmat((0:sx + 2*extra_x - 1), [N 1]);
    for i = 1:N
        I1(r(i),c(i),:) = g * reshape(I(YI(i,:),XI(i,:),:), [R2,b]);
%         I1(r(i),c(i),:) = sum(g1 .* reshape(I(YI(i,:),XI(i,:),:), [R2,b]), 1);
    end
else % use parallel toolbox to speed up, it seems to make it slower
    I1 = zeros(N,b);
%     generate_data = @(tmp) I;
%     const_I = WorkerObjWrapper(generate_data);
    const_I = parallel.pool.Constant(I);
    parfor i = 1:N
        [r,c] = ind2sub([rows,cols],i);
        yi = (r-1)*sy + 1 : r*sy + 2*extra_y;
        xi = (c-1)*sx + 1 : c*sx + 2*extra_x;
        rect = const_I.Value(yi,xi,:);
        I1(i,:) = sum(g1 .* reshape(rect, [R2,b]), 1);
    end
    I1 = reshape(I1,[rows,cols,b]);
end

% g = g(:);
% C = construct_C(rows,cols,[s2,s1]);
% 
% C1 = cat(2,C{:});
% G = g2G(C1,g);
% I2 = reshape(I1, [size(I1,1)*size(I1,2), size(I1,3)]);
% I3 = G * I2;
% I1 = reshape(I3, [rows, cols, size(I3,2)]);
