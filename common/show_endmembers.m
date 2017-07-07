function show_endmembers(M2,wl,names,options)
if nargin < 3
    M = size(M2,1);
    names = cell(1,M);
    for i = 1:M
        names{i} = ['endmember ',num2str(i)];
    end
end
if nargin < 4
    options = [];
end

bbl = parse_param(options,'bbl',[]);
figure;
if ~iscell(M2)
    lines = show_spectra_in_one_plot(wl, M2, bbl);
    set_line_styles(lines);
    legend(names);
else
    plot_row = ceil(sqrt(length(M2)));
    plot_col = ceil(length(M2)/plot_row);
    for i = 1:length(M2)
        subplot(plot_row,plot_col,i);
        show_spectra_in_one_plot(wl, M2{i}, bbl);
        xlabel({'Wavelength (micrometer)',names{i}});
    end
end
end

function lines = show_spectra_in_one_plot(wl, M2, bbl)
M2(:,~bbl) = NaN;
lines = plot(repmat(wl(:),1,size(M2,1)),M2','LineWidth',1);
ylim([0 1]);
xlabel('Wavelength (micrometer)');
ylabel('Reflectance');
end

function set_line_styles(lines)
M = length(lines);
colors = distinguishable_colors(M);
for i = 1:M
    set(lines(i),'color',colors(i,:));
end

end