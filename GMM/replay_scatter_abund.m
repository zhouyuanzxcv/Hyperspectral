function replay_scatter_abund(frames_scatter, frames_abund, fps)
%REPLAY_SCATTER_ABUND Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    fps = 5;
end

concat = 0;
if length(frames_scatter) == length(frames_abund)
    scatter1 = frames_scatter(1);
    abund1 = frames_abund(1);
    if isempty(scatter1.colormap) && isempty(abund1.colormap)
        concat = 1;
    end
end

if ~concat
    scatter_fh = figure('name','Replay evolution of the Gaussians');
    movie(scatter_fh,frames_scatter,1,fps);
    abund_fh = figure('name','Replay evolution of the abundances');
    movie(abund_fh,frames_abund,1,fps);
    return;
end

N = length(frames_scatter);
disp(['Replay ',num2str(N), ' frames']);
for i = 1:N
    scatter_data = frames_scatter(i).cdata;
    abund_data = frames_abund(i).cdata;
    new_data = cat(2,scatter_data,abund_data);
    new_frames(i).cdata = new_data;
    new_frames(i).colormap = [];
end

fh = figure('name','Replay evolution of the Gaussians and the abundances');
pos = get(fh,'position');
pos(1) = pos(1) - pos(3)/2;
pos(3) = pos(3) * 2;
set(fh,'position',pos);
movie(fh,new_frames,1,fps);