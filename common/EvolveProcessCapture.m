classdef EvolveProcessCapture
    % Shows the evolution process of spectral unmixing algorithms.
    properties
      ScatterOpts
      AbundOpts
      FramesScatter
      FramesAbund
    end
    methods
        function obj = EvolveProcessCapture()
        end
        
        function obj = begin(obj,max_iter)
            scatter_opts = [];
            scatter_opts.hs = [];
            scatter_opts.update = 0;

            abund_opts = [];
            abund_opts.imhs = [];
            abund_opts.update = 0;
            
            frames_scatter(max_iter) = struct('cdata',[],'colormap',[]);
            frames_abund(max_iter) = struct('cdata',[],'colormap',[]);
            
            obj.ScatterOpts = scatter_opts;
            obj.AbundOpts = abund_opts;
            obj.FramesScatter = frames_scatter;
            obj.FramesAbund = frames_abund;
        end
        
        function obj = end(obj,end_iter)
            obj.FramesScatter = obj.FramesScatter(1:end_iter);
            obj.FramesAbund = obj.FramesAbund(1:end_iter);
        end
        
        function obj = updateEndmVar(obj,Y,A,rows,cols,w_jk,mu_jk,sigma_jk,names,iter)
            obj = updateScatterPlotEndmVar(obj,Y,w_jk,mu_jk,sigma_jk,names,iter);
            obj = updateAbundanceMaps(obj,A,rows,cols,iter);

            drawnow
        end
        
        function obj = updateEndmFixed(obj,Y,A,rows,cols,R,names,iter)
            obj = updateScatterPlotEndmFixed(obj,Y,R,names,iter);
            obj = updateAbundanceMaps(obj,A,rows,cols,iter);
            
            drawnow
        end
        
        function obj = updateScatterPlotEndmVar(obj,Y,w_jk,mu_jk,sigma_jk,names,iter)
            scatter_opts = obj.ScatterOpts;
            frames_scatter = obj.FramesScatter;
            
            [~,hs,axes_h,scatter_fh] = endmember_scatter_plot_end_var( ...
                Y,w_jk,mu_jk,sigma_jk,names,scatter_opts);
            title(axes_h, ['Scatter plot of iteration ',num2str(iter)]);
            frames_scatter(iter) = getframe(scatter_fh);

            if ~scatter_opts.update
                pos = get(scatter_fh,'Position');
                pos(1) = pos(1) - pos(3) / 2 - 1;
                set(scatter_fh,'Position',pos);
            end
            scatter_opts.update = 1;
            scatter_opts.hs = hs;
            scatter_opts.axes_h = axes_h;
            scatter_opts.fig_h = scatter_fh;
            
            obj.ScatterOpts = scatter_opts;
            obj.FramesScatter = frames_scatter;
        end
        
        function obj = updateScatterPlotEndmFixed(obj,Y,R,names,iter)
            scatter_opts = obj.ScatterOpts;
            frames_scatter = obj.FramesScatter;
            
            [~,hs,axes_h,scatter_fh] = endmember_scatter_plot( ...
                Y,R,names,scatter_opts);
            title(axes_h, ['Scatter plot of iteration ',num2str(iter)]);
            frames_scatter(iter) = getframe(scatter_fh);

            if ~scatter_opts.update
                pos = get(scatter_fh,'Position');
                pos(1) = pos(1) - pos(3) / 2 - 1;
                set(scatter_fh,'Position',pos);
            end
            scatter_opts.update = 1;
            scatter_opts.hs = hs;
            scatter_opts.axes_h = axes_h;
            scatter_opts.fig_h = scatter_fh;
            
            obj.ScatterOpts = scatter_opts;
            obj.FramesScatter = frames_scatter;
        end

        
        function obj = updateAbundanceMaps(obj,A,rows,cols,iter)
            abund_opts = obj.AbundOpts;
            frames_abund = obj.FramesAbund;
            
            [imhs,abund_fh] = show_abundances(A,rows,cols,'',0,0,abund_opts);
            set(abund_fh, 'name', ['Abundances of iteration ',num2str(iter)]);
            frames_abund(iter) = getframe(abund_fh);

            if ~abund_opts.update
                pos = get(abund_fh,'Position');
                pos(1) = pos(1) + pos(3) / 2 + 1;
                set(abund_fh,'Position',pos);
            end
            abund_opts.update = 1;
            abund_opts.imhs = imhs;
            abund_opts.fig_h = abund_fh;
            
            obj.AbundOpts = abund_opts;
            obj.FramesAbund = frames_abund;
        end
    end
end