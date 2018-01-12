function add_path(onedrive, dropbox)
if nargin < 2
    if ispc
        onedrive = 'D:/OneDrive/';
        dropbox = 'D:/Dropbox/';
    elseif ismac
        onedrive = '/Users/yuanzhou/OneDrive/';
        dropbox = '/Users/yuanzhou/Dropbox/';
    elseif isunix
        dropbox = '/cise/homes/yuan/Dropbox/';
        onedrive = '/cise/homes/yuan/OneDrive/';
    end
end

if ispc
    addpath('C:/Program Files/Xpdf/bin64');
end

addpath([dropbox,'YuanHyperspectral/code/Grendel/Graphcut']);

addpath([onedrive,'Code/Matlab/Utility/distinguishable_colors']);
addpath([onedrive,'Code/Matlab/Utility/export_fig']);
addpath([onedrive,'Code/Matlab/Utility/gridLegend']);
addpath([onedrive,'Code/Matlab/Utility/subtightplot']);
addpath([onedrive,'Code/Matlab/Utility/textprogressbar']);
addpath([onedrive,'Code/Matlab/Utility']);

s = version;
if str2double(s(1)) <= 7
    addpath([onedrive,'Code/Matlab/Utility/compatibility']);
end

addpath([dropbox,'YuanHyperspectral/code/common']);
addpath([dropbox,'YuanHyperspectral/code/common/load_image']);
