function show_uncertainty_range(M1, M2_1, var_amts, var_dirs, wl, names, options)
if nargin < 7
    options = [];
end

[M,B] = size(M2_1);
row = parse_param(options,'plot_row',floor(sqrt(M)));
col = ceil(M/row);

use_uncertainty = ~isempty(var_amts);

linewidth = 1.5;

for j = 1:size(M2_1,1)
    if all(isnan(M2_1(j,:)))
        M2_1(j,:) = -1;
    end
end

figure('name','Compare two endmember spectra');
for j = 1:M
    subplot(row,col,j); hold on;
    if ~isempty(M1)
        h_gt = plot(wl,M1(j,:),'k-','LineWidth',linewidth);
    end
    h_uncer = -1;
    if use_uncertainty
        h_uncer = plot(wl,2*var_amts(j)*var_dirs(j,:) + M2_1(j,:),...
            ':','LineWidth',linewidth,'Color',[0 0.5 0]);
        plot(wl,-2*var_amts(j)*var_dirs(j,:) + M2_1(j,:),...
            ':','LineWidth',linewidth,'Color',[0 0.5 0]);
    end
    h_algo = plot(wl,M2_1(j,:),'r--','LineWidth',linewidth);

    ylim([0 1]);

    xlabel({'Wavelength (micrometer)',['(',char(96+j),') ',names{j}]},'fontsize',12);
    ylabel('Reflectance','fontsize',12);
end

if h_uncer == -1;
    if ~isempty(M1)
        h_leg = legend([h_gt,h_algo],'Ground Truth','Algorithm');
    else
        h_leg = legend(h_algo,'Algorithm');
    end
else
    if ~isempty(M1)
        h_leg = legend([h_gt,h_algo,h_uncer],'Ground Truth','Algorithm','Uncertainty Range');
    else
        h_leg = legend([h_algo,h_uncer],'Algorithm','Uncertainty Range');
    end
end
set(h_leg,'Location','Best','Orientation','horizontal','fontsize',12);

