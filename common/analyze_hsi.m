function varargout = analyze_hsi(varargin)
% ANALYZE_HSI MATLAB code for analyze_hsi.fig
%      ANALYZE_HSI, by itself, creates a new ANALYZE_HSI or raises the existing
%      singleton*.
%
%      H = ANALYZE_HSI returns the handle to a new ANALYZE_HSI or the handle to
%      the existing singleton*.
%
%      ANALYZE_HSI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYZE_HSI.M with the given input arguments.
%
%      ANALYZE_HSI('Property','Value',...) creates a new ANALYZE_HSI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyze_hsi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyze_hsi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyze_hsi

% Last Modified by GUIDE v2.5 24-Jan-2016 15:26:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyze_hsi_OpeningFcn, ...
                   'gui_OutputFcn',  @analyze_hsi_OutputFcn, ...
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


% --- Executes just before analyze_hsi is made visible.
function analyze_hsi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyze_hsi (see VARARGIN)

% Choose default command line output for analyze_hsi
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analyze_hsi wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Doesn't help
% set(handles.uitoggletool_probe,'state','off');

dcm_obj = datacursormode;
set(dcm_obj,'UpdateFcn',{@my_data_cursor_updatefcn,hObject})

function txt = my_data_cursor_updatefcn(~,event_obj,hObject)
% Customizes text of data tips
handles = guidata(hObject);
pos = get(event_obj,'Position');
ind = get(event_obj, 'DataIndex');
txt = num2str(pos);

% show spectra
spectra = squeeze(handles.I(pos(2),pos(1),:));
wl = handles.wl';
plot(handles.axes_spectra,wl,spectra);
ylim(handles.axes_spectra,[0 1]);

% show the other cursor
handles = delete_ui_handles(handles,'other_cursor');
hold(handles.axes_gt,'on');
other_cursor = plot(handles.axes_gt,pos(1),pos(2),'s','MarkerSize',5,...
    'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',2);
hold(handles.axes_gt,'off');
handles.other_cursor = other_cursor;

% update status
set(handles.text_status, 'string', get_gt_name_by_position(handles,pos));

guidata(hObject, handles);

function gt_name = get_gt_name_by_position(handles, pos)
gt_vec = squeeze(handles.A_gt(pos(2),pos(1),:));
gt_name = strjoin(handles.gt_names(gt_vec==1),', ');

% --- Outputs from this function are returned to the command line.
function varargout = analyze_hsi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function pavia_Callback(hObject, eventdata, handles)
% hObject    handle to pavia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[I,R_gt,A_gt,names,wl,rgb] = load_pavia_university([],0);
tick_labels = {'Asphalt','Meadows','Gravel','Trees','Metal',...
    'Soil','Bitumen','Bricks','Shadows'};
M = size(A_gt,3);
colors = distinguishable_colors(M);
I_gt = A_gt2I_gt(A_gt);

handles = load_image(rgb,I,R_gt,A_gt,names,wl,I_gt,colors,tick_labels,handles);
guidata(hObject,handles);


function handles = load_image(rgb,I,R_gt,A_gt,names,wl,I_gt,colors,tick_labels,handles)
% show rgb image and the ground truth
imshow(rgb,'parent',handles.axes_image);
show_gt_image(I_gt,colors,handles.axes_gt,tick_labels);
show_gt_spectra(R_gt,wl,colors,handles.axes_gt_spectra,tick_labels)

% load gt names into the listbox
display_names = [{'All'}, names];
set(handles.listbox_gt_names,'String',display_names,...
	'Value',1);
set(handles.popupmenu_change_gt,'String',names,'Value',1);

% save handles
handles.rgb = rgb;
handles.I = I;
handles.R_gt = R_gt;
handles.A_gt = A_gt;
handles.I_gt = I_gt;
handles.tick_labels = tick_labels;
handles.gt_colors = colors;
handles.gt_names = names;
handles.wl = wl;

function show_gt_image(I_gt,colors,axes_gt,tick_labels)
I_gt = uint8(I_gt);
colors = [1 1 1;colors];
M = length(tick_labels);
im = findobj(axes_gt,'type','image');
if isempty(im) || ~isequal(size(get(im,'CData')), size(I_gt))
    imshow(I_gt,colors,'parent',axes_gt);
else
    set(im,'CData',I_gt);
    colormap(colors);
end


function show_gt_spectra(R_gt,wl,colors,axes_gt_spectra,tick_labels)
M = length(tick_labels);
cla(axes_gt_spectra);
hold(axes_gt_spectra,'on');
for i = 1:M
    plot(axes_gt_spectra,wl,R_gt(i,:),'LineWidth',0.5,'Color',colors(i,:));
end
hold(axes_gt_spectra,'off');
ylim(axes_gt_spectra,[0 1]);

% xlabel('Wavelength (micrometer)','fontsize',5);
% ylabel('Reflectance','fontsize',5);
lgh = legend(axes_gt_spectra, tick_labels);
set(lgh,'fontsize',5);


function edit_crop_pos_Callback(hObject, eventdata, handles)
% hObject    handle to edit_crop_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_crop_pos as text
%        str2double(get(hObject,'String')) returns contents of edit_crop_pos as a double


% --- Executes during object creation, after setting all properties.
function edit_crop_pos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_crop_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_crop.
function pushbutton_crop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
crop_pos = handles.crop_pos;
set(handles.edit_crop_pos,'string','x1 y1 x2 y2');
handles = delete_ui_handles(handles,'roi_rectangle');
handles.crop_pos = [];

x1 = crop_pos(1);
y1 = crop_pos(2);
x2 = crop_pos(3);
y2 = crop_pos(4);

I = handles.I;
A_gt = handles.A_gt;
rgb = handles.rgb;
gt_names1 = handles.gt_names;
I_gt = handles.I_gt;
gt_colors1 = handles.gt_colors;
tick_labels1 = handles.tick_labels;

I1 = I(y1:y2,x1:x2,:);
A_gt1 = A_gt(y1:y2,x1:x2,:);
rgb1 = rgb(y1:y2,x1:x2,:);
M = size(A_gt,3);

[R_gt1,index_to_be_removed] = get_gt_mean_spectra(I1,A_gt1,1:M);

A_gt1(:,:,index_to_be_removed) = [];
R_gt1(index_to_be_removed,:) = [];
gt_names1(index_to_be_removed) = [];
tick_labels1(index_to_be_removed) = [];
gt_colors1(index_to_be_removed,:) = [];

% Reset the I_gt values to range from 1 to M
I_gt1 = A_gt2I_gt(A_gt1);

handles = load_image(rgb1,I1,R_gt1,A_gt1,gt_names1,handles.wl,I_gt1, ...
    gt_colors1,tick_labels1,handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function select_roi_Callback(hObject, eventdata, handles)
% hObject    handle to select_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fh = figure;imshow(handles.rgb);
[roi,x1,y1,x2,y2] = selectroi(getimage);
close(fh);
handles.crop_pos = [x1,y1,x2,y2];
handles = update_crop_display(handles);
guidata(hObject,handles);


% --- Executes on key press with focus on edit_crop_pos and none of its controls.
function edit_crop_pos_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_crop_pos (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    handles.crop_pos = str2num(get(hObject,'string'));
    handles = update_crop_display(handles);
    guidata(hObject,handles);
end

function handles = update_crop_display(handles)
roi = handles.crop_pos;
set(handles.edit_crop_pos,'string',num2str(handles.crop_pos));

handles = delete_ui_handles(handles,'roi_rectangle');
handles.roi_rectangle = draw_roi_rectangle(handles.axes_image, roi);

function rect_h = draw_roi_rectangle(axes_h, roi)
x = [roi(1),roi(3),roi(3),roi(1)];
y = [roi(2),roi(2),roi(4),roi(4)];
x = [x,x(1)];
y = [y,y(1)];

hold(axes_h,'on');
rect_h = plot(axes_h,x,y,'r-','Linewidth',2);
hold(axes_h,'off');

function handles = delete_ui_handles(handles, ui_str)
if isfield(handles,ui_str)
    try
        delete(handles.(ui_str));
    catch me
        disp(['Error to delete ui handles ',num2str(handles.(ui_str))]);
    end
    handles.(ui_str) = [];
end


% --- Executes on selection change in listbox_gt_names.
function listbox_gt_names_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_gt_names (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_gt_names contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_gt_names

if strcmp(get(handles.figure1,'SelectionType'),'normal') % If single left click
    index_selected = get(handles.listbox_gt_names,'Value');
    I_gt = handles.I_gt;
    colors = handles.gt_colors;
    M = size(colors,1);
    tick_labels = handles.tick_labels;
    axes_gt = handles.axes_gt;
    if index_selected == 1
        show_gt_image(I_gt,colors,axes_gt,tick_labels);
    else
        ind = index_selected - 1;
        for j = 1:M
            if j ~= ind
                I_gt(I_gt==j) = 0;
            end
        end
        show_gt_image(I_gt,colors,axes_gt,tick_labels);
    end
    handles.I_gt_current = I_gt;
elseif strcmp(get(hObject,'SelectionType'),'open') % If double click
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function listbox_gt_names_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_gt_names (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit_crop_pos,'string','x1 y1 x2 y2');
handles = delete_ui_handles(handles,'roi_rectangle');

guidata(hObject,handles);


% --- Executes on selection change in popupmenu_modify_size.
function popupmenu_modify_size_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_modify_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_modify_size contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_modify_size


% --- Executes during object creation, after setting all properties.
function popupmenu_modify_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_modify_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_add_gt.
function radiobutton_add_gt_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_add_gt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_add_gt


% --- Executes on button press in radiobutton_remove_gt.
function radiobutton_remove_gt_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_remove_gt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_remove_gt


% --------------------------------------------------------------------
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rgb = handles.rgb;
I = handles.I;
R_gt = handles.R_gt;
A_gt = handles.A_gt;
names = handles.gt_names;
wl = handles.wl;
I_gt = handles.I_gt;
gt_colors = handles.gt_colors;
tick_labels = handles.tick_labels;

[savefile, pathname] = uiputfile('*.mat','Save Dataset As');

set(gcf,'Pointer','watch');
drawnow;
save(fullfile(pathname, savefile),'rgb','I','R_gt','A_gt','names', ...
    'wl','I_gt','gt_colors','tick_labels');
set(gcf,'Pointer','arrow');

% --------------------------------------------------------------------
function menu_dataset_Callback(hObject, eventdata, handles)
% hObject    handle to menu_dataset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.mat','Select the dataset file');
if isequal(filename,0)
    return;
end

load(fullfile(pathname,filename));
handles = load_image(rgb,I,R_gt,A_gt,names,wl,I_gt,gt_colors,tick_labels,handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function uitoggletool_probe_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_probe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disable_other_togglebuttons_in_toolbar(handles,hObject);


function disable_other_togglebuttons_in_toolbar(handles,current)
all_togglebuttons = [handles.uitoggletool_probe, ...
    handles.uitoggletool_zoom_in, handles.uitoggletool_zoom_out, ...
    handles.uitoggletool_pan, handles.uitoggletool_modify_gt];
% This does not actually disables the actions
for togglebutton = all_togglebuttons
    if togglebutton ~= current
%         set(togglebutton,'state','off');
    end
end


% --------------------------------------------------------------------
function uitoggletool_zoom_in_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_zoom_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disable_other_togglebuttons_in_toolbar(handles,hObject);


% --------------------------------------------------------------------
function uitoggletool_zoom_out_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_zoom_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disable_other_togglebuttons_in_toolbar(handles,hObject);


% --------------------------------------------------------------------
function uitoggletool_pan_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disable_other_togglebuttons_in_toolbar(handles,hObject);


% --------------------------------------------------------------------
function uitoggletool_modify_gt_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_modify_gt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disable_other_togglebuttons_in_toolbar(handles,hObject);

% handles.old_button_down_fcn = get(gcf,'WindowButtonDownFcn');
% handles.old_button_motion_fcn = get(gcf,'WindowButtonMotionFcn');
% handles.old_button_up_fcn = get(gcf,'WindowButtonUpFcn');

set(gcf,'WindowButtonDownFcn',@my_button_down_fcn);
set(gcf,'WindowButtonMotionFcn',@my_button_motion_fcn);
set(gcf,'WindowButtonUpFcn',@my_button_up_fcn);

% set(gcf,'Pointer','crosshair');


% --------------------------------------------------------------------
function uitoggletool_modify_gt_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_modify_gt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if isfield(handles,'old_button_down_fcn') && ...
%         isfield(handles,'old_button_motion_fcn')
set(handles.figure1,'WindowButtonDownFcn',[]);
set(handles.figure1,'WindowButtonMotionFcn',[]);
set(handles.figure1,'WindowButtonUpFcn',[]);
% end
% set(gcf,'Pointer','arrow');
handles = delete_ui_handles(handles, 'modify_circle');
handles = delete_ui_handles(handles, 'reference_circle');
handles = delete_ui_handles(handles, 'prepare_circle');

axes_hs = [handles.axes_image, handles.axes_gt];
for axes_h = axes_hs
    circles = findobj(axes_h,'type','rectangle');
    delete(circles);
end

guidata(hObject,handles);

function my_button_down_fcn(src,eventdata)
% get(gcf,'SelectionType')
handles = guidata(src);
if ~test_hit_in_image(handles) || get(handles.listbox_gt_names,'Value') == 1
    return;
end
handles.mouse_left_pressed = 1;
handles = delete_ui_handles(handles, 'prepare_circle');
handles = perform_modification(handles);
guidata(src,handles);

function my_button_motion_fcn(src,eventdata)
handles = guidata(src);
if ~test_hit_in_image(handles)
    return;
end

curPt = get(handles.axes_image, 'CurrentPoint');
x = curPt(1,1);
y = curPt(1,2);
pos = [round(x),round(y)];

% show spectra
spectra = squeeze(handles.I(pos(2),pos(1),:));
wl = handles.wl';
plot(handles.axes_spectra,wl,spectra);
ylim(handles.axes_spectra,[0 1]);
% update status
curr_gt = get_gt_name_by_position(handles, pos);
curr_pos = ['(',num2str(x,'% 10.1f'),', ',num2str(y,'% 10.1f'),')'];
set(handles.text_status, 'string', [curr_gt, ' ',curr_pos]);

if isfield(handles,'mouse_left_pressed') && handles.mouse_left_pressed && ...
        get(handles.listbox_gt_names,'Value') > 1;
    handles = perform_modification(handles);
else
    handles = draw_prepare_circle(handles);
end
guidata(src,handles);

function in = test_hit_in_image(handles)
curPt = get(handles.axes_image, 'CurrentPoint');
% Get changed pixels
x = round(curPt(1,1));
y = round(curPt(1,2));
if x >= 1 && x <= size(handles.I_gt,2) && y >= 1 && y <= size(handles.I_gt,1)
    in = 1;
else
    in = 0;
end

function handles = perform_modification(handles)
curPt = get(handles.axes_image, 'CurrentPoint');
% Get changed pixels
x = curPt(1,1);
y = curPt(1,2);
radius = get_pen_radius(handles.popupmenu_modify_size);
[X,Y] = meshgrid([round(x-radius):round(x+radius)], ...
    [round(y-radius):round(y+radius)]);
X = [X(:),Y(:)];
d = sqrt(sum((X - repmat([x,y],size(X,1),1)).^2, 2));
X1 = X(d <= radius, :);
X1(X1(:,1) < 1 | X1(:,1) > size(handles.I_gt,2), :) = [];
X1(X1(:,2) < 1 | X1(:,2) > size(handles.I_gt,1), :) = [];
I_gt = handles.I_gt_current;
ind = sub2ind(size(I_gt),X1(:,2),X1(:,1));

% draw new gt
gt_ind = get(handles.listbox_gt_names,'Value') - 1;

colors = handles.gt_colors;
tick_labels = handles.tick_labels;

if gt_ind == 0
    errordlg('Select one kind of ground truth in the listbox to modify pixels');
end

if ~isempty(ind) && gt_ind > 0
    if get(handles.radiobutton_add_gt,'Value')
        I_gt(ind) = gt_ind;
        show_gt_image(I_gt,colors,handles.axes_gt,tick_labels);
    elseif get(handles.radiobutton_remove_gt,'Value')
        I_gt(ind) = 0;
        show_gt_image(I_gt,colors,handles.axes_gt,tick_labels);
    elseif get(handles.radiobutton_change_gt,'Value')
        change_to_ind = get(handles.popupmenu_change_gt,'Value');
        if  change_to_ind ~= gt_ind
            selected_pixels = I_gt(ind);
            ind = ind(selected_pixels ~= 0);
            I_gt(ind) = 0;
            A_gt1 = handles.A_gt(:,:,change_to_ind);
            A_gt1(ind) = 1;
            handles.A_gt(:,:,change_to_ind) = A_gt1;
            show_gt_image(I_gt,colors,handles.axes_gt,tick_labels);
        end
    end
    handles.gt_modified = 1;
end
handles.I_gt_current = I_gt;

handles = draw_modify_circle(handles);

function radius = get_pen_radius(radius_h)
radius_ind = get(radius_h,'Value'); 
radius_strs = get(radius_h, 'String'); 
radius = str2num(radius_strs{radius_ind});

function handles = draw_modify_circle(handles)
curPt = get(handles.axes_image, 'CurrentPoint');
x = curPt(1,1);
y = curPt(1,2);

radius = get_pen_radius(handles.popupmenu_modify_size);

handles = delete_ui_handles(handles, 'modify_circle');
handles.modify_circle = draw_circle(handles.axes_image,x,y,radius,'r');

handles = delete_ui_handles(handles, 'reference_circle');
handles.reference_circle = draw_circle(handles.axes_gt,x,y,radius,'r');


function handles = draw_prepare_circle(handles)
curPt = get(handles.axes_image, 'CurrentPoint');
x = curPt(1,1);
y = curPt(1,2);

radius = get_pen_radius(handles.popupmenu_modify_size);

handles = delete_ui_handles(handles, 'prepare_circle');
handles.prepare_circle = draw_circle(handles.axes_image,x,y,radius,'g');

handles = delete_ui_handles(handles, 'reference_circle');
handles.reference_circle = draw_circle(handles.axes_gt,x,y,radius,'g');


function circle_h = draw_circle(axes_h,x,y,radius,color)
hold(axes_h,'on');
circle_h = rectangle('Position', [x-radius,y-radius,2*radius,2*radius], ...
    'Curvature',[1,1],'FaceColor',color,'Parent',axes_h);
hold(axes_h,'off');


function my_button_up_fcn(src,eventdata)
handles = guidata(src);
if ~handles.mouse_left_pressed || ~handles.gt_modified || ~test_hit_in_image(handles)
    return;
end

set(gcf,'Pointer','watch');
drawnow;

handles.mouse_left_pressed = 0;
handles.gt_modified = 0;
handles = delete_ui_handles(handles, 'modify_circle');
handles = draw_prepare_circle(handles);

handles = update_gt(handles);

set(gcf,'Pointer','arrow');
guidata(src,handles);

function handles = update_gt(handles)
I = handles.I;
wl = handles.wl;
I_gt = handles.I_gt;
A_gt = handles.A_gt;
colors = handles.gt_colors;
tick_labels = handles.tick_labels;
M = size(A_gt,3);
gt_ind = get(handles.listbox_gt_names,'Value') - 1;

I_gt_current1 = (handles.I_gt_current ~= 0);
A_gt(:,:,gt_ind) = double(I_gt_current1);
for j = 1:M
    if j ~= gt_ind
        A_gt1 = A_gt(:,:,j);
        A_gt1(I_gt_current1) = 0;
        A_gt(:,:,j) = A_gt1;
    end
end

I_gt = A_gt2I_gt(A_gt);

[R_gt,index_to_be_removed] = get_gt_mean_spectra(I,A_gt);
show_gt_spectra(R_gt,wl,colors,handles.axes_gt_spectra,tick_labels);
handles.A_gt = A_gt;
handles.R_gt = R_gt;
handles.I_gt = I_gt;


% --------------------------------------------------------------------
function uitoggletool_probe_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool_probe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = delete_ui_handles(handles,'other_cursor');
guidata(hObject,handles);


% --- Executes on selection change in popupmenu_change_gt.
function popupmenu_change_gt_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_change_gt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_change_gt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_change_gt


% --- Executes during object creation, after setting all properties.
function popupmenu_change_gt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_change_gt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gt_name_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gt_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gt_name as text
%        str2double(get(hObject,'String')) returns contents of edit_gt_name as a double


% --- Executes during object creation, after setting all properties.
function edit_gt_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gt_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_add_gt_name.
function pushbutton_add_gt_name_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_gt_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gt_name = get(handles.edit_gt_name,'string');
if isempty(gt_name)
    errordlg('Can not add an empty ground truth name');
    return;
end

idx = get(handles.listbox_gt_names,'value');
selected_str = get_selected_string_from_listbox(handles.listbox_gt_names);

R_gt = handles.R_gt;
A_gt = handles.A_gt;
gt_names = handles.gt_names;
I_gt = handles.I_gt;
gt_colors = handles.gt_colors;
tick_labels = handles.tick_labels;

% Modify the background data
gt_names = insert_by_dim(gt_names, {gt_name}, idx, 2);
R_gt = insert_by_dim(R_gt, zeros(1,size(R_gt,2)), idx, 1);
A_gt = insert_by_dim(A_gt, zeros(size(A_gt,1),size(A_gt,2)), idx, 3);
I_gt = A_gt2I_gt(A_gt);
tick_labels = insert_by_dim(tick_labels, {gt_name}, idx, 2);
M = size(A_gt,3);

colors = distinguishable_colors(M);
for j = 1:M
    candidate_color = colors(j,:);
    find_matched_color = 0;
    for k = 1:size(gt_colors,1)
        if all(candidate_color == gt_colors(k,:))
            find_matched_color = 1;
            break;
        end
    end
    if ~find_matched_color
        break;
    end
end
gt_colors = insert_by_dim(gt_colors, candidate_color, idx, 1);
        
handles.R_gt = R_gt;
handles.A_gt = A_gt;
handles.I_gt = I_gt;
handles.gt_colors = gt_colors;
handles.tick_labels = tick_labels;
handles.gt_names = gt_names;

% Change UI display
listbox_strs = [{'All'} gt_names];
for j = 1:length(listbox_strs)
    if strcmp(listbox_strs{j},selected_str)
        break;
    end
end
refresh_gt_names_display(handles, gt_names, j);

guidata(hObject,handles);

function refresh_gt_names_display(handles, gt_names, idx)
set(handles.listbox_gt_names,'value',idx);
listbox_strs = [{'All'} gt_names];
set(handles.listbox_gt_names,'string',listbox_strs);
set(handles.popupmenu_change_gt,'string',gt_names);
show_gt_spectra(handles.R_gt, handles.wl, handles.gt_colors, ...
    handles.axes_gt_spectra, handles.tick_labels);

function selected_str = get_selected_string_from_listbox(listbox)
idx = get(listbox,'value');
listbox_strs = get(listbox,'string');
selected_str = listbox_strs{idx};

% --- Executes on button press in pushbutton_remove_gt_name.
function pushbutton_remove_gt_name_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_gt_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selected_str = get_selected_string_from_listbox(handles.listbox_gt_names);
idx = get(handles.listbox_gt_names,'value');
idx = idx - 1;

if strcmp(selected_str,'All') || idx == 0
    errordlg('Must select one specific ground truth to delete.');
    return;
end

% Modify background data
R_gt = handles.R_gt;
A_gt = handles.A_gt;
gt_names = handles.gt_names;
I_gt = handles.I_gt;
gt_colors = handles.gt_colors;
tick_labels = handles.tick_labels;

R_gt(idx,:) = [];
A_gt(:,:,idx) = [];
gt_names(idx) = [];
gt_colors(idx,:) = [];
tick_labels(idx) = [];
I_gt = A_gt2I_gt(A_gt);

handles.R_gt = R_gt;
handles.A_gt = A_gt;
handles.I_gt = I_gt;
handles.gt_colors = gt_colors;
handles.tick_labels = tick_labels;
handles.gt_names = gt_names;

% refresh display
refresh_gt_names_display(handles, gt_names, idx);

guidata(hObject,handles);

listbox_gt_names_Callback(handles.listbox_gt_names, [], handles);


% --- Executes on button press in pushbutton_change_gt_name.
function pushbutton_change_gt_name_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_change_gt_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gt_name = get(handles.edit_gt_name,'string');
if isempty(gt_name)
    errordlg('Can not change the selected ground truth to an empty name.');
    return;
end

idx = get(handles.listbox_gt_names,'value');
idx = idx - 1;
if idx == 0
    errordlg('Must select one specific ground truth to change name.');
    return;
end

% Modify background data
R_gt = handles.R_gt;
gt_colors = handles.gt_colors;
gt_names = handles.gt_names;
tick_labels = handles.tick_labels;

gt_names{idx} = gt_name;
tick_labels{idx} = gt_name;

handles.gt_names = gt_names;
handles.tick_labels = tick_labels;

% refresh display
refresh_gt_names_display(handles, gt_names, idx+1);

guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
imshow(handles.rgb);

if isfield(handles,'crop_pos') && ~isempty(handles.crop_pos);
    draw_roi_rectangle(gca, handles.crop_pos);
end

% --------------------------------------------------------------------
function menu_gt_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_gt_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
show_gt_image(handles.I_gt, handles.gt_colors, gca, handles.tick_labels);
M = length(handles.tick_labels);
% show colorbar
cbh = colorbar;
set(cbh,'Location','West');
% pos = get(cbh,'Position');
% pos(1) = 0.9;
% set(cbh,'Location','manual');
% set(cbh,'Position',pos);
set(cbh,'YTickLabel',handles.tick_labels);
set(cbh,'YDir','reverse');

s = version;
if str2double(s(1)) <= 7
    set(cbh,'YTick',(1:M)');
    set(cbh,'YLim',[0.5 M+0.5]);
else
    set(cbh,'ticks',(0.5:1:M+0.5)');
    set(cbh,'YLim',[0 M]);
end

% --------------------------------------------------------------------
function menu_gt_spectra_Callback(hObject, eventdata, handles)
% hObject    handle to menu_gt_spectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
subplot(2,2,1);
show_gt_spectra(handles.R_gt, handles.wl, handles.gt_colors, ...
    gca, handles.tick_labels);
lines = findobj(gca,'type','line');
for line = lines
    set(line,'linewidth',2);
end

xlabel('Wavelength (micrometer)','fontsize',12);
ylabel('Reflectance','fontsize',12);
lgh = legend(gca, handles.gt_names);
set(lgh,'fontsize',12);


% --------------------------------------------------------------------
function menu_gt_Callback(hObject, eventdata, handles)
% hObject    handle to menu_gt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_gt_hist_Callback(hObject, eventdata, handles)
% hObject    handle to menu_gt_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Y,A,rows,cols] = reshape_hsi(handles.I, handles.A_gt);
hist_end_var(Y, A, handles.gt_names);

% --------------------------------------------------------------------
function menu_reset_colors_Callback(hObject, eventdata, handles)
% hObject    handle to menu_reset_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
M = length(handles.gt_names);
colors2 = distinguishable_colors(M);

rgb = im2double(handles.rgb);
if 1
    [colors1,index_to_be_removed] = get_gt_mean_spectra(rgb, handles.A_gt);
    [P,colors] = permute_endmembers(colors1, colors2);
elseif 0
    [colors,index_to_be_removed] = get_gt_mean_spectra(rgb, handles.A_gt);
else
    colors = colors2;
end

handles.gt_colors = colors;

show_gt_spectra(handles.R_gt, handles.wl, colors, ...
    handles.axes_gt_spectra, handles.tick_labels);

guidata(hObject,handles);

listbox_gt_names_Callback(handles.listbox_gt_names, [], handles);


% --------------------------------------------------------------------
function menu_script_Callback(hObject, eventdata, handles)
% hObject    handle to menu_script (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.m','Select the load script file');
if isequal(filename,0)
    return;
end
filepath = fullfile(pathname,filename);
[pathstr,name,ext] = fileparts(filepath);
% current_dir = pwd;
% cd(pathstr);
eval(['[I,R_gt,A_gt,names,wl,rgb] = ',name,';']);
% cd(current_dir);

% If no ground truth, set default one
if isempty(R_gt) && isempty(A_gt)
    % set default number of endmembers
    default_num_endmembers = 1;
    
    [Y,A_gt,rows,cols] = reshape_hsi(I,A_gt);
    [idx,c] = kmeans(Y,default_num_endmembers);
    R_gt = c;
    A_gt = I_gt2A_gt(uint8(reshape(idx, rows, cols)));
    for i = 1:size(R_gt,1)
        names{i} = ['endmember',num2str(i)];
    end
end

% set other missing information
if ~exist('rgb','var') || isempty(rgb)
    rgb = retrieve_rgb(I,wl);
end
I_gt = A_gt2I_gt(A_gt);
M = size(A_gt,3);
gt_colors = distinguishable_colors(M);
tick_labels = names;

% load image
handles = load_image(rgb,I,R_gt,A_gt,names,wl,I_gt,gt_colors,tick_labels,handles);
guidata(hObject,handles);
