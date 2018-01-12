function [mapping, hs, axes_h, fig_h] = endmember_scatter_plot_impl(Y,M,names,options)
%ENDMEMBER_SCATTER_PLOT_IMPL Summary of this function goes here
%   Detailed explanation goes here
Xs = [];
update = 0;
hs = cell(1, M + 1);
axes_h = [];
fig_h = [];

dimension_num = 2;
colors = distinguishable_colors(M);

fcn_draw_endmembers = [];

if nargin > 3 && isstruct(options)
    arg_set = fieldnames(options);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=options.',arg_set{i},';']);
    end
end

% B = size(mus{1},2);
options.dimension_num = dimension_num;
use_specific_projection = parse_param(options, 'use_specific_projection', []);
if isempty(use_specific_projection)
    [mappedX, mapping] = pca(Y,dimension_num);
else
    mapping = options.use_specific_projection;
    mapping.M = mapping.M(:,1:dimension_num);
    mapping.lambda = mapping.lambda(1:dimension_num);
    mappedX = (Y - repmat(mapping.mean, size(Y,1), 1)) * mapping.M;
end

%% Show the pixel cloud
if ~update
    hs = [];
    fig_h = figure;
    if dimension_num == 2
        hs{1} = scatter(mappedX(:,1), mappedX(:,2), 10, '.', ...
            'MarkerEdgeColor', [0.6 0.6 0.6]); 
    elseif dimension_num == 3
        hs{1} = scatter3(mappedX(:,1), mappedX(:,2), mappedX(:,3), 10, '.', ...
            'MarkerEdgeColor', [0.6 0.6 0.6]); 
    end
    hold on; axis equal;
    axes_h = gca;
end

%% Show the specified pixels in another color
if ~isempty(Xs)
    for i = 1:length(Xs)
        marker_color = colors(i,:) / 2;
        if norm(marker_color - colors(i,:), 1) < 0.1
            colors(i,:) = colors(i,:) + 0.5;
        end
%         marker_color = min(ones(1,3), colors(i,:) + [0.5 0.5 0.5]);
        mappedX = (Xs{i} - repmat(mapping.mean, size(Xs{i},1), 1)) * mapping.M;
        if dimension_num == 2
            scatter(axes_h, mappedX(:,1), mappedX(:,2), 5, '+', ...
                'MarkerEdgeColor', marker_color);
        elseif dimension_num == 3
            scatter3(axes_h, mappedX(:,1), mappedX(:,2), mappedX(:,3), 5, '+', ...
                'MarkerEdgeColor', marker_color);
        end
    end
end

if update
    for j = 2:length(hs)
        delete(hs{j});
        hs{j} = [];
    end
end

hs = fcn_draw_endmembers(axes_h, mapping, colors, hs, options);

% set_marker(hs);
names = {'Pixels',names{:}};

hs_legend = zeros(1,length(hs));
for j = 1:length(hs)
    hs_legend(j) = hs{j}(1);
end
legend(hs_legend, names,'Location','best');

function set_marker(hs)
set(hs(1),'Marker','s');
if length(hs) > 1
    set(hs(2),'Marker','*');
end
if length(hs) > 2
    set(hs(3),'Marker','x');
end
if length(hs) > 3
    set(hs(4),'Marker','o');
end


