function [ output_args ] = hist_wavelength_reflectance(Y,A_gt,names,wl,options)
%HIST_WAVELENGTH_REFLECTANCE Summary of this function goes here
%   Detailed explanation goes here
if nargin < 5
    options = [];
end

[N,B] = size(Y);
M = size(A_gt,2);

thresh = parse_param(options, 'thresh', 0.99);
reflectance_bin_num = parse_param(options, 'reflectance_bin_num', 51);
w_jk = parse_param(options, 'w_jk', []);
mu_jk = parse_param(options, 'mu_jk', []);
sigma_jk = parse_param(options, 'sigma_jk', []);
rows = parse_param(options, 'plot_row', ceil(sqrt(M)));
cols = ceil(M/rows);

% sort the components by prior probability
if ~isempty(w_jk)
    for j = 1:length(w_jk)
        [sorted_w,ind] = sort(w_jk{j},'descend');
        sorted_mu = mu_jk{j}(ind,:);
        sorted_sigma = sigma_jk{j}(:,:,ind);
        w_jk{j} = sorted_w;
        mu_jk{j} = sorted_mu;
        sigma_jk{j} = sorted_sigma;
    end
end

reflectance_bins = linspace(0,1,reflectance_bin_num);

X_coord = repmat(wl,[reflectance_bin_num,1]);
Y_coord = (reflectance_bins([2:end,end]) + reflectance_bins)' / 2;
Y_coord = repmat(Y_coord,[1,length(wl)]);

colors = distinguishable_colors(7);
colors([1,4],:) = [];
colors([3,4],:) = colors([4,3],:);

figure,

for j = 1:M
    subplot(rows,cols,j);
    % show the background histogram image
    Z = zeros(reflectance_bin_num,B);
    Y1 = Y(A_gt(:,j) >= thresh,:);
    if ~isempty(Y1)
        for k = 1:B
            bincounts = histc(Y1(:,k), reflectance_bins);
            Z(:,k) = bincounts;
        end
        normalized_Z = Z / sum(Z(:));
        im = normalized_Z * B;
        im = uint8(round(im*255));
        cm = colormap('gray');
        cm = flipud(cm);
        imagesc(X_coord(1,:)',Y_coord(:,1),im); axis xy;
        colormap(cm);
    end

    xlabel({'Wavelength (micrometer)',names{j}});
    ylabel('Reflectance');
    
    % add the gmm curves
    line_hs = [];
    line_str = cell(0);
    if ~isempty(w_jk)
        for k = 1:length(w_jk{j})
            line_str{end+1} = num2str(w_jk{j}(k),2);
            center = mu_jk{j}(k,:)';
            cov = sigma_jk{j}(:,:,k);
            cov = (1/2) * (cov+cov'); % force symmetry
            [V,D] = eig(cov);
            d = sqrt(diag(D));
            [sd,ind] = max(d);
            dir = V(:,ind);
            hold on;
            center_h = plot_gaussian_curve(wl',center,'--',colors(k,:));
            plot_gaussian_curve(wl',center + 2*dir*sd,':',colors(k,:));
            plot_gaussian_curve(wl',center - 2*dir*sd,':',colors(k,:));
            line_hs(end+1) = center_h;
        end
        legend(line_hs,line_str,'location','best');
    end
    set(gca,'xlim',[min(wl),max(wl)],'ylim',[0 1]);
end
        
function h = plot_gaussian_curve(wl,y,linestyle,color)
h = plot(wl, y, 'linestyle', linestyle, 'color', color,'linewidth',1);


