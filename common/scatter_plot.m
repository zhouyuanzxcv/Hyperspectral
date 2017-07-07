function [mapping, hs, axes_h, fig_h] = scatter_plot(Y,names,options)
%ENDMEMBER_SCATTER_PLOT Summary of this function goes here
%   Detailed explanation goes here
M = length(names);
Xs = [];
colors = distinguishable_colors(M);
dimension_num = 2;

if nargin > 2 && isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end

[mappedX, mapping] = pca(Y,dimension_num);

fig_h = figure;
hs = cell(0);
if dimension_num == 2
    scatter_h = scatter(mappedX(:,1), mappedX(:,2), 10, '.', ...
        'MarkerEdgeColor', [0.6 0.6 0.6]);
elseif dimension_num == 3
    scatter_h = scatter3(mappedX(:,1), mappedX(:,2), mappedX(:,3), 10, '.', ...
        'MarkerEdgeColor', [0.6 0.6 0.6]);
end

hs(end+1) = {scatter_h};
hold on; axis equal;
axes_h = gca;

for i = 1:length(Xs)
    marker_color = colors(i,:);
    mappedX = (Xs{i} - repmat(mapping.mean, size(Xs{i},1), 1)) * mapping.M;
    if dimension_num == 2
        scatter_h = scatter(axes_h, mappedX(:,1), mappedX(:,2), 5, '+', ...
            'MarkerEdgeColor', marker_color);
    elseif dimension_num == 3
        scatter_h = scatter3(axes_h, mappedX(:,1), mappedX(:,2), ...
            mappedX(:,3), 5, '+', 'MarkerEdgeColor', marker_color);
    end
    hs(end+1) = {scatter_h};
end

% show legend
names = {'Pixels',names{:}};

hs_legend = zeros(1,length(hs));
for j = 1:length(hs)
    hs_legend(j) = hs{j}(1);
end
legend(hs_legend, names,'Location','best');
