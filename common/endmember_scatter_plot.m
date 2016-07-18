function [mapping, hs, axes_h, fig_h] = endmember_scatter_plot(Y,R,names,options)
%ENDMEMBER_SCATTER_PLOT Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    options = [];
end

[M,B,p] = size(R);

options.R = R;
options.fcn_draw_endmembers = @draw_endmembers;

[mapping, hs, axes_h, fig_h] = endmember_scatter_plot_impl(Y,M,names,options);


function hs = draw_endmembers(axes_h, mapping, colors, hs, options)
R = options.R;
if isempty(R)
    hs = [];
    return;
end

[M,B,p] = size(R);

for i = 1:p
    proj_R(:,:,i) = (R(:,:,i) - repmat(mapping.mean,M,1))*mapping.M;
end

% for i = 1:M
%     pt_labels{i} = num2str(i);
% end

% The LineWidth will determine the edge thickness of the markers
for j = 1:M
    hs1 = plot(axes_h,squeeze(proj_R(j,1,:)),squeeze(proj_R(j,2,:)),'s',...
        'MarkerSize',8);
    set(hs1,'color',colors(j,:),'linewidth',1.5);
    hs{j+1} = hs1;
end
% labelpoints(proj_R(:,1,1),proj_R(:,2,1), pt_labels);
