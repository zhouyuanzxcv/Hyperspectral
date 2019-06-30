function add_path(dropbox)
if nargin < 1
    if ispc
        dropbox = 'C:/Dropbox/';
    elseif ismac
        dropbox = '/Users/yuanzhou/Dropbox/';
    elseif isunix
        dropbox = '/cise/homes/yuan/Dropbox/';
    end
end

if ispc
    addpath('C:/Program Files/Xpdf/bin64');
end

addpath([dropbox,'YuanHyperspectral/code/Grendel/Graphcut']);

addpath([dropbox,'Code/Matlab/Utility/distinguishable_colors']);
addpath([dropbox,'Code/Matlab/Utility/export_fig']);
addpath([dropbox,'Code/Matlab/Utility/gridLegend']);
addpath([dropbox,'Code/Matlab/Utility/subtightplot']);
addpath([dropbox,'Code/Matlab/Utility/textprogressbar']);
addpath([dropbox,'Code/Matlab/Utility']);

s = version;
if str2double(s(1)) <= 7
    addpath([dropbox,'Code/Matlab/Utility/compatibility']);
end

addpath([dropbox,'YuanHyperspectral/code/common']);
addpath([dropbox,'YuanHyperspectral/code/common/load_image']);
