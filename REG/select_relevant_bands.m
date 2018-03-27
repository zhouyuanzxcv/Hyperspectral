function [I_sel, bands_sel] = select_relevant_bands(I,wl,mode)
%SELECT_RELAVENT_BANDS Summary of this function goes here
%   Detailed explanation goes here

% we assume wl has unit micrometer
if ischar(mode)
    if strcmp(mode, 'full')
        range = [min(wl),max(wl)];
    elseif strcmp(mode, 'multispectral')
        range = [0.4,0.9];
    elseif strcmp(mode, 'color')
        range = [0.4,0.8];
    elseif strcmp(mode, 'panchromatic')
        range = [0.4,1]; % 600nm
    end
else
    range = mode;
end

bands_sel = find(wl >= range(1) & wl <= range(2));
I_sel = I(:,:,bands_sel);

% wavelength_width = mean(wl(2:end) - wl(1:end-1));
% 
% size_srf = round((range / wavelength_width) / 2);
% 
% sel_mat = zeros(length(wl),sel_mat_cols);
% 
% for k = 1:size(sel_mat,2)
%     start_ind = rgb_ind(k)-size_srf;
%     end_ind = rgb_ind(k)+size_srf;
%     start_ind = max(start_ind, 1);
%     end_ind = min(end_ind, size(sel_mat,1));
%     sel_mat(start_ind:end_ind, k) = 1;
% end
% 
% bands_sel = find(sum(sel_mat,2) > 0);
% I_sel = I(:,:,bands_sel);


end

