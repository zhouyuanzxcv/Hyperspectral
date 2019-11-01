function varargout = hyper_rgb_comparison(varargin)
% HYPER_RGB_COMPARISON MATLAB code for hyper_rgb_comparison.fig
%      HYPER_RGB_COMPARISON, by itself, creates a new HYPER_RGB_COMPARISON or raises the existing
%      singleton*.
%
%      H = HYPER_RGB_COMPARISON returns the handle to a new HYPER_RGB_COMPARISON or the handle to
%      the existing singleton*.
%
%      HYPER_RGB_COMPARISON('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HYPER_RGB_COMPARISON.M with the given input arguments.
%
%      HYPER_RGB_COMPARISON('Property','Value',...) creates a new HYPER_RGB_COMPARISON or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hyper_rgb_comparison_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hyper_rgb_comparison_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hyper_rgb_comparison

% Last Modified by GUIDE v2.5 02-Mar-2017 14:23:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hyper_rgb_comparison_OpeningFcn, ...
                   'gui_OutputFcn',  @hyper_rgb_comparison_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before hyper_rgb_comparison is made visible.
function hyper_rgb_comparison_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hyper_rgb_comparison (see VARARGIN)

% Choose default command line output for hyper_rgb_comparison
addpath('reg_fft');
addpath('transform');
addpath('optimization');

handles.output = hObject;

%% load images
rho = 0;
% dataset = 'salton_sea';
dataset = get(hObject,'UserData');

switch dataset
    case 'pavia'
        U2 = zeros(80,50);
        V2 = zeros(80,50);
        bbl = true(103,1);
        rho = 3;
        load('reg_pavia_dataset.mat');
        load('reg_result_pavia.mat');
        I = I1;
    case 'salton_sea'
%         rho = ceil((20/16.9) * mean(s2) / 2);

        load('salton_sea_roi.mat');
%         rgb1 = imread('salton_sea_color.png');
        rgb1 = imread('salton_sea_roi_3.2015.jpg');

        load('reg_result_salton_sea.mat');
%         load('reg_result_salton_sea_rigid.mat');
%         load('reg_result_salton_sea_MI.mat');
%         load('reg_result_salton_sea_MI_bspline.mat');
        
        bbl = logical(bbl);
    case 'neon_sjer'
        load('neon_sjer_roi.mat');
        rgb1 = imread('2013_SJER_AOP_Camera_sample.tif');
        load('reg_result_neon_sjer.mat');
        bbl = logical(bbl);
    otherwise
end

% show hyperspectral image
imshow(uint8(retrieve_rgb(I, wl) * 255),'Parent',handles.axes_hyper);
show_reg_warp(handles.axes_hyper, zeros(size(U2)), zeros(size(V2)));

% show hyperspectral image warping field
quiver(handles.axes_vector_field,U2,V2); 
set(handles.axes_vector_field,'YDir','reverse');

% show warped hyperspectral image
options = [];
options.isRigid = 0;
options.useTranslationField = 1;
options.U = U2;
options.V = V2;
options.rho = rho;

I(:,:,~bbl) = 0;
I_hyper_reg = transform(I,eye(3),[1,1],1e-3,size(I,2),size(I,1),options);
I_hyper_reg(:,:,~bbl) = NaN;
I(:,:,~bbl) = NaN;
imshow(uint8(retrieve_rgb(I_hyper_reg, wl) * 255),'Parent',handles.axes_warped_hyper);
show_reg_warp(handles.axes_warped_hyper, U2, V2);

% show_reg_warp(retrieve_rgb(I, wl));
% show_reg_warp(retrieve_rgb(I_hyper_reg, wl), U2, V2);

% transform color image
[rows,cols,B] = size(I);
s = s2;
T = T2;
sigma = sigma2;

s1 = ceil(s(1)); % scale_x
s2 = ceil(s(2)); % scale_y

options.isRigid = 1;
[I1,rgb,g,extra] = transform(double(rgb1),T,s,sigma,cols,rows,options);

% show transformed color image
imshow(uint8(I1),'Parent',handles.axes_transformed_color_image);

% show residual image and SRF
[residual,H1,wl_sel] = calc_residual_image(I_hyper_reg(:,:,bbl), I1, wl(bbl));
imshow(uint8(residual*2),'Parent',handles.axes_residual);
plot(handles.axes_srf, wl_sel, H1(2:end,1), 'r-', wl_sel,H1(2:end,2), ...
    'g-', wl_sel, H1(2:end,3), 'b-');

% show rgb image
imshow(uint8(rgb),'Parent',handles.axes_rgb);
top_left_pt = extra+0.5;
bottom_right_pt = size(rgb);
bottom_right_pt = bottom_right_pt(2:-1:1) - extra + 0.5;
draw_roi_rectangle(handles.axes_rgb,[top_left_pt,bottom_right_pt],'r--');

% show PSF
mesh(handles.axes_psf, g);

handles.s = [s1,s2];
handles.I = I;
handles.I_hyper_reg = I_hyper_reg;
handles.rgb = rgb;
handles.rgb_transformed = I1;
handles.psf = g;
handles.wl = wl;
handles.extra = extra;
handles.I_view_spectra = I_hyper_reg; % switched over popup menu
handles.U = U2;
handles.V = V2;
handles.H1 = H1;
handles.residual = residual;
set(handles.popupmenu_spectra_selector, 'Value', 1);

%% Add data cursor action
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes hyper_rgb_comparison wait for user response (see UIRESUME)
% uiwait(handles.figure_hyper_rgb);
dcm_obj = datacursormode(handles.figure_hyper_rgb);
set(dcm_obj,'UpdateFcn',{@my_data_cursor_updatefcn,hObject})


function [residual,H1,wl_sel] = calc_residual_image(I, I1, wl)
[I_sel,bands_sel] = select_relevant_bands(I,wl,'color');
wl_sel = wl(bands_sel)';

Y1 = reshape_hsi(I_sel);
Y1 = [ones(size(Y1,1),1),Y1];

X = reshape_hsi(I1);
% options = [];
% options.Y1 = Y1;
% options.Y1Y = Y1'*Y1;
% options.Y1Yinv = inv(options.Y1Y);
options = [];
options.lambda = 1e-3;

H1 = solve_for_H1(X, Y1, options);

Y2 = Y1*H1;
residual = sqrt(mean((X - Y2).^2, 2));
residual = reshape(residual, [size(I,1) size(I,2)]);

% calculate other metric e.g. correlation coefficients (corrcoef)
X_m = mean(X,1);
Y_m = mean(Y2,1);
X3 = X - repmat(X_m,[size(X,1),1]);
Y3 = Y2 - repmat(Y_m,[size(X,1),1]);
CC = sum(sum(X3.*Y3, 2)) / sqrt(sum(sum(X3.^2, 2)) * sum(sum(Y3.^2, 2)));
disp('Correlation coefficient is ')
CC

function txt = my_data_cursor_updatefcn(~,event_obj,hObject)
%% Customizes text of data tips
handles = guidata(hObject);
pos = get(event_obj,'Position');
ind = get(event_obj, 'DataIndex');
image = get(event_obj,'Target');
data = get(image,'CData');
pixel = data(pos(2),pos(1),:);
txt = ['Pos: ',num2str(pos),', Data: ',num2str(pixel)];

%% show spectra
spectra = squeeze(handles.I_view_spectra(pos(2),pos(1),:));
wl = handles.wl';
plot(handles.axes_spectra,wl,spectra,'linewidth',1);
ylim(handles.axes_spectra,[0 1]);
% xlabel('Micrometer');
% ylabel('Reflectance');
handles.axes_spectra.UserData = [];
handles.axes_spectra.UserData.wl = wl;
handles.axes_spectra.UserData.spectra = spectra;

%% show the other cursors
handles = delete_ui_handles(handles,'other_cursor');
handles = delete_ui_handles(handles,'psf_cursor');

s = handles.s;

[rect_exact, rect_psf] = calculate_correspondence_region(pos,s, ...
    handles.extra, handles.rgb, handles.I);

other_cursor = draw_roi_rectangle(handles.axes_rgb,rect_exact,'r-');
psf_cursor = draw_psf_circle(handles.axes_rgb,rect_psf,'g-');

handles.cursor_pos = pos;
handles.other_cursor = other_cursor;
handles.psf_cursor = psf_cursor;

guidata(hObject, handles);


function handles = delete_ui_handles(handles, ui_str)
if isfield(handles,ui_str)
    try
        delete(handles.(ui_str));
    catch me
        disp(['Error to delete ui handles ',num2str(handles.(ui_str))]);
    end
    handles.(ui_str) = [];
end


% --- Outputs from this function are returned to the command line.
function varargout = hyper_rgb_comparison_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_spectra_selector.
function popupmenu_spectra_selector_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_spectra_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_spectra_selector contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_spectra_selector
value = get(handles.popupmenu_spectra_selector, 'Value');
if value == 1
    handles.I_view_spectra = handles.I_hyper_reg;
elseif value == 2
    handles.I_view_spectra = handles.I;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_spectra_selector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_spectra_selector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function menu_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuitem_spectra_Callback(hObject, eventdata, handles)
% hObject    handle to menuitem_spectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
wl = handles.axes_spectra.UserData.wl;
spectra = handles.axes_spectra.UserData.spectra;
fh = figure; 
plot(wl,spectra,'linewidth',1);
ylim([0,1]);
xlabel('Wavelength (micrometer)');
ylabel('Reflectance');
set(fh,'pos', [488.2, 555.4, 267.2, 187.2]);


% --------------------------------------------------------------------
function menuitem_nonrigid_translation_Callback(hObject, eventdata, handles)
% hObject    handle to menuitem_nonrigid_translation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% show hyperspectral image warping field
U2 = handles.U;
V2 = handles.V;
quiver_distortion_field(U2, V2);


% --------------------------------------------------------------------
function menuitem_profile_envi_Callback(hObject, eventdata, handles)
% hObject    handle to menuitem_profile_envi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = handles.cursor_pos; % center point
wl = handles.wl;

styles = {'r-','g-','b-','k-'};
[rgb_roi,rects,spectras] = get_profile_envi_data(pos,handles);

% show results
gap = 0.01;
marg_h = 0.15;
marg_w = 0.15;

figure;
subtightplot(2,1,1,gap,marg_h,marg_w);
imshow(uint8(rgb_roi));
for i = 1:length(rects)
    draw_psf_circle(gca, rects{i}, styles{i});
end

subtightplot(2,1,2,gap,marg_h,marg_w);
hold on;
for i = 1:length(spectras)
    plot(wl,spectras{i},styles{i},'linewidth',1);
end
xlabel('Wavelength (micrometer)');
ylabel('Reflectance');
ylim([0 1]);


function [rgb_roi,rects,spectras] = get_profile_envi_data(pos,handles)
pos1 = pos;
pos1(1) = pos1(1) - 1; % left point
pos2 = pos;
pos2(1) = pos2(1) + 1; % right point
pos3 = pos;
pos3(1) = pos3(1) + 2; % next right point

poses = {pos1,pos,pos2,pos3};

s = handles.s;
extra = handles.extra;
rgb = handles.rgb;
I = handles.I;

rects = cell(1,length(poses));
spectras = cell(1,length(poses));

for i = 1:length(poses)
    pos = poses{i};
    [rect_exact, rect_psf] = calculate_correspondence_region(pos,s,extra,rgb,I);
    rects{i} = rect_psf;
end

rect_roi_extra = [6*s(1),1*s(2)];
rect_roi = [floor(rects{1}(1:2) - rect_roi_extra), ...
    ceil(rects{end}(3:4) + rect_roi_extra)];
rgb_roi = rgb(rect_roi(2):rect_roi(4),rect_roi(1):rect_roi(3),:);
roi_shift = rect_roi(1:2) - 1;
roi_shift = repmat(roi_shift, [1 2]);

for i = 1:length(poses)
    rects{i} = rects{i} - roi_shift;
end

I_hyper_reg = handles.I_hyper_reg;
for i = 1:length(poses)
    pos = poses{i};
    spectra = squeeze(I_hyper_reg(pos(2),pos(1),:));
    spectras{i} = spectra;
end


% --------------------------------------------------------------------
function menuitem_all_images_Callback(hObject, eventdata, handles)
% hObject    handle to menuitem_all_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.I;
wl = handles.wl;
I_hyper_reg = handles.I_hyper_reg;
rgb_transformed = handles.rgb_transformed;
rgb = handles.rgb;

fh = figure;

subtightplot(2,2,1); % orginal hyperspectral image
imshow(uint8(retrieve_rgb(I, wl) * 255));

subtightplot(2,2,2); % registered color image
imshow(uint8(rgb));

subtightplot(2,2,3); % warped hyperspectral image
imshow(uint8(retrieve_rgb(I_hyper_reg, wl) * 255));

subtightplot(2,2,4); % transformed color image
imshow(uint8(rgb_transformed));

set(fh,'position',[488 342 464 420]);


% --------------------------------------------------------------------
function menuitem_residual_image_Callback(hObject, eventdata, handles)
% hObject    handle to menuitem_residual_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
residual = handles.residual;
figure,imshow(uint8(residual*2));
disp('The mean of residual map is ');
mean(residual(:))
set(gcf,'position',[398,95,923,594]);

% --------------------------------------------------------------------
function menuitem_custom_Callback(hObject, eventdata, handles)
% hObject    handle to menuitem_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
locations = [9,3;15,10;22,45;30,25];
wl = handles.wl;

% show results
figure;

gap = [0.01,0.05]; % (vertical,horizontal) gap between axes
marg_h = 0.2;
marg_w = 0.1;
for j = 1:size(locations,1)
    pos = locations(j,:);
    styles = {'r-','g-','b-','k-'};
    [rgb_roi,rects,spectras] = get_profile_envi_data(pos,handles);
    [rgb_roi,rects] = pad_rgb_roi(rgb_roi, rects);

    subtightplot(2,4,j,gap,marg_h,marg_w);
    imshow(uint8(rgb_roi));
    for i = 1:length(rects)
        draw_psf_circle(gca, rects{i}, styles{i});
    end

    subtightplot(2,4,j+4,gap,marg_h,marg_w);
    hold on;
    for i = 1:length(spectras)
        plot(wl,spectras{i},styles{i},'linewidth',1);
    end
    xlabel('Wavelength (micrometer)');
    ylabel('Reflectance');
    ylim([0 1]);
end

set(gcf,'position',[325,268,1078,386]);

function [rgb_roi1,rects1] = pad_rgb_roi(rgb_roi, rects)
[rows,cols,b] = size(rgb_roi);
rgb_roi1 = ones(rows*3,cols,b)*255;
rgb_roi1(2*rows+1:end,:,:) = rgb_roi;
rects1 = rects;
for i = 1:length(rects1)
    rect = rects1{i};
    rect(2) = rect(2) + 2*rows;
    rect(4) = rect(4) + 2*rows;
    rects1{i} = rect;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% utility functions %%%%%%%%%%%%%%%%%%%%%%%%%
function [rect_exact, rect_psf] = calculate_correspondence_region(pos,s,extra,rgb,I)
% check the size is correct
assert(extra(1)*2 + size(I,2)*s(1) == size(rgb,2));
assert(extra(2)*2 + size(I,1)*s(2) == size(rgb,1));

% imshow has 1,2,3,... as coordinates for the centers of the pixels. Hence
% [10.5, 20.5] will encompass the pixels we are interested in (11,12,...,20)
% for s = 10, pos = 2.
center = [(pos(1)-0.5)*s(1), (pos(2)-0.5)*s(2)] + 0.5;
top_left_pt = [center(1)-s(1)/2, center(2)-s(2)/2];
bottom_right_pt = [center(1)+s(1)/2, center(2)+s(2)/2];

% shift by the margin size
top_left_pt = top_left_pt + extra;
bottom_right_pt = bottom_right_pt + extra;
rect_exact = [top_left_pt,bottom_right_pt];

% expand by margin to indicate the range covered by the PSF
top_left_pt = top_left_pt - extra;
bottom_right_pt = bottom_right_pt + extra;
rect_psf = [top_left_pt,bottom_right_pt];


function quiver_distortion_field(U2, V2)
[rows,cols] = size(U2);

x = (0.5:cols-0.5);
y = (0.5:rows-0.5);
[X,Y] = meshgrid(x,y);
figure,quiver(X,Y,U2,V2);
set(gca,'YDir','reverse');
set(gca,'xlim',[0 cols]);
set(gca,'ylim',[0 rows]);


function rect_h = draw_roi_rectangle(axes_h, roi, style)
% roi = (x1,y1,x2,y2) with (x1,y1) being the top left corner
x = [roi(1),roi(3),roi(3),roi(1)];
y = [roi(2),roi(2),roi(4),roi(4)];
x = [x,x(1)];
y = [y,y(1)];

hold(axes_h,'on');
rect_h = plot(axes_h,x,y,style,'Linewidth',1);
hold(axes_h,'off');


function circle_h = draw_psf_circle(axes_h, roi, style)
px = roi(1);
py = roi(2);
dx = roi(3) - roi(1);
dy = roi(4) - roi(2);

hold(axes_h,'on');
circle_h = rectangle('Position',[px py dx dy],'Curvature',[1,1], ...
    'edgecolor',style(1),'linestyle',style(2),'parent',axes_h,'linewidth',1);
hold(axes_h,'off');

