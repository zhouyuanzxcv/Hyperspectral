function [mapping, hs, axes_h, fig_h] = endmember_scatter_plot_end_var(Y,ws,mus,sigmas,names,options)
%ENDMEMBER_SCATTER_PLOT Summary of this function goes here
%   Detailed explanation goes here
M = length(ws);
options.ws = ws;
options.mus = mus;
options.sigmas = sigmas;
options.fcn_draw_endmembers = @draw_endmembers;

[mapping, hs, axes_h, fig_h] = endmember_scatter_plot_impl(Y,M,names,options);


function hs = draw_endmembers(axes_h, mapping, colors, hs, options)
ws = options.ws;
mus = options.mus;
sigmas = options.sigmas;
M = length(ws);

if isempty(ws{1})
    return;
end

[proj_R, proj_Cov] = project_gmm(mus, sigmas, mapping);


for j = 1:M
    h1 = [];
    for k = 1:size(proj_R{j},1)
        h1(k) = plot_gaussian_ellipsoid(proj_R{j}(k,:), ...
            proj_Cov{j}(:,:,k), 2, 50, axes_h);
        set(h1(k),'color',colors(j,:),'linewidth',2);     
    end
    hs{j+1} = h1;
%     labelpoints(proj_R{j}(1,1), proj_R{j}(1,2), pt_labels{j});
end