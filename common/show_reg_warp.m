function show_reg_warp(I,U,V,grid_size)
%SHOW_REG_WARP Summary of this function goes here
%   Detailed explanation goes here
% Input:
%   I - image matrix or handle to the axes
%   (U,V) - translation field v(x) = (U(x),V(x)), s.t. T(x) = x + v(x)
%   grid_size - number of points per dimension

if nargin < 1
    I = imread('pout.tif');
end

if nargin < 3
    [rows,cols,b] = size(I);
    U = zeros(rows,cols);
    V = zeros(rows,cols);
end

if nargin < 4
    grid_size = 20;
end

[rows,cols] = size(U);


x1 = linspace(1, cols, grid_size);
y1 = linspace(1, rows, grid_size);
[X,Y] = meshgrid(x1,y1);
U1 = interp2((1:cols),(1:rows),U,X,Y);
V1 = interp2((1:cols),(1:rows),V,X,Y);
X_target = X - U1; % X_target + U1 = X;
Y_target = Y - V1; % Y_target + V1 = Y;
x = X_target(:);
y = Y_target(:);

[W,M] = image2graph(ones(length(y1),length(x1)),0.05,1e-9);

% [A,J2,A2] = create_adjacency_graph(X,'epsball',1.01,1);

if ishandle(I)
    ax_h = I;
elseif ~isempty(I)
    figure,imshow(I);
    ax_h = gca;
end
hold(ax_h, 'on');

for i = 1:length(x)
    neighbor_inds = find(M(i,:));
    for j = 1:length(neighbor_inds)
        plot(ax_h,[x(i),x(neighbor_inds(j))],[y(i),y(neighbor_inds(j))],'-b',...
            'linewidth', 1);
    end
end

end

