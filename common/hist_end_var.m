function hist_end_var(Y,A,names,thresh,options)
%HIST_END_VAR Display histogram of pure pixels.
% Input:
%   Y - N by B pixel data
%   A - N by M abundance maps
%   names - ground truth names

if nargin < 4
    thresh = 0.99;
end

if nargin < 5
    options = [];
end

show_approx = parse_param(options, 'show_approx', 0);
w_jk = parse_param(options, 'w_jk', []);
mu_jk = parse_param(options, 'mu_jk', []);
sigma_jk = parse_param(options, 'sigma_jk', []);
legend_names = parse_param(options, 'legend_names', []);
line_styles = parse_param(options, 'line_styles', {'--','-.',':'});
colors = parse_param(options, 'colors',[1 0 0;0.5 0.5 0;0 0.75 0.75]);

[N,B] = size(Y);
[~,M] = size(A);
num_bins = 20;

% [mappedX, mapping] = pca(Y,2);

figure,

hs_legend = [];
rows = ceil(sqrt(M));
cols = ceil(M/rows);
for j = 1:M
    subplot(rows,cols,j);
    Y1 = Y(A(:,j) >= thresh,:);
    if isempty(Y1)
        disp(['Warning: The set of samples for endmember ',num2str(j),' is empty.']);
        xlabel(names{j});
        continue;
    end
    [X1, mapping] = pca(Y1,1);
    
%     num_bins = round(size(X1,1) / 10);
%     num_bins = min(max(num_bins,10),30);
    [bar_h, bin_values] = display_hist(X1,num_bins);
%     hs_legend(1) = bar_h;
    hs_legend = [];

    xlabel(names{j});
    
    if show_approx
        for k = 1:size(mu_jk,1)
            [proj_mu, proj_sigma] = project_gmm(mu_jk(k,:), sigma_jk(k,:), mapping);
            value = calc_log_gmm(bin_values, w_jk{k,j}, proj_mu{j}, proj_sigma{j});
            y = (bin_values(2) - bin_values(1)) * exp(value);
            hold on; 
            
            lh = plot(bin_values, y, 'LineStyle', line_styles{k}, ...
                'color', colors(k,:),  'LineWidth', 2);
            hs_legend = [hs_legend, lh];
        end
    end
end

if ~isempty(legend_names)
    legend(hs_legend,legend_names,'Orientation','horizontal','Location','best');
end

function [bar_h, binValues] = display_hist(x, num_bins)
[counts, binValues] = hist(x, num_bins);
normalizedCounts = counts / sum(counts);

bar_h = bar(binValues, normalizedCounts, 'barwidth', 1);
% xlabel('Input Value');
ylabel('Probability');

binValues = binValues';