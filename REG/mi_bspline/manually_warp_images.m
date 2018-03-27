function varargout = manually_warp_images(varargin)
% MANUALLY_WARP_IMAGES, is a GUI which shows the controlpoint grid used
% during registration of two images, allowing the user to drag control 
% points to get a better initial alignement before pressing the start 
% registration button.
%
% Example,
%   Istatic=im2double(imread('images/prostate1.png'));
%   Imoving=im2double(imread('images/prostate2.png'));
%
%   manually_warp_images(Imoving,Istatic);
%
% note: on line 410 of this code you can find the options used 
% for the registration...
%
% Function is written by D.Kroon University of Twente (January 2010)


% Edit the above text to modify the response to help manually_warp_images

% Last Modified by GUIDE v2.5 07-Jan-2010 14:37:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manually_warp_images_OpeningFcn, ...
                   'gui_OutputFcn',  @manually_warp_images_OutputFcn, ...
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


% --- Executes just before manually_warp_images is made visible.
function manually_warp_images_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manually_warp_images (see VARARGIN)

% Choose default command line output for manually_warp_images
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes manually_warp_images wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if (isempty(varargin))
    error('manually_warp_images:inputs','no images set');
else
    data.Istatic=varargin{1};
    data.Imoving=varargin{2};
end

try
    functionname='manually_warp_images.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions'])
    addpath([functiondir '/functions_affine'])
    addpath([functiondir '/functions_nonrigid'])
catch me
    disp(me.message);
end

% Disable warning
warning('off', 'MATLAB:maxNumCompThreads:Deprecated')

[data.Istatic,data.Imoving]=images2samesize(data.Istatic,data.Imoving);

% Calculate max refinements steps
MaxItt=min(floor(log2([size(data.Imoving,1) size(data.Imoving,2)]/4)));

% set b-spline grid spacing in x and y direction
data.Spacing=[2^MaxItt 2^MaxItt];
data.Uniform=make_init_grid(data.Spacing,[size(data.Imoving,1) size(data.Imoving,2)]);
data.Grid=data.Uniform;
data.axesgrid_handle=handles.axes1;
data.axespreview1_handle=handles.axes2;
data.axespreview2_handle=handles.axes3;
data.pointselectedhandle=[];
data.padding=50;
data.Similarity='mi';
data.viewoption=0;
[data.Iclass,data.Imin,data.Imax,data.Imoving,data.Istatic]=images2double(data.Imoving,data.Istatic);
data.Itotal=zeros([size(data.Imoving,1)+data.padding*2 size(data.Imoving,2)+data.padding*2 3]);
data.imshowhandle=imshow(data.Itotal,'Parent',data.axesgrid_handle);
data.imshowpreview1handle=imshow(data.Imoving,'Parent',data.axespreview1_handle);
data.imshowpreview2handle=imshow(data.Istatic,'Parent',data.axespreview2_handle);
data.gridpointhandle=[];
data.handles=handles;
if(size(data.Istatic,3)==1)
    data.Jstatic=zeros([size(data.Istatic) 3]);
    data.Jstatic(:,:,2)=0.5*data.Istatic; data.Jstatic(:,:,3)=data.Istatic; 
else
    data.Jstatic=data.Istatic;
end

data.GUItimer = timer('TimerFcn', 'manually_warp_images(''TimerFcn'',gcf,[],guidata(gcf));','Period',1,'ExecutionMode','fixedDelay'); 
setMyData(data);
start(data.GUItimer);
bsplinetranformimage();
set_menu_checks();
showslices();
showgrid();

function [Istatic,Imoving]=images2samesize(Istatic,Imoving)
% Resize the moving image to fit the static image
if(sum(size(Istatic)-size(Imoving))~=0)
    Imoving = imresize(Imoving,[size(Istatic,1) size(Istatic,2)],'bicubic');
end


function [Iclass,Imin,Imax,Imoving,Istatic]=images2double(Imoving,Istatic)
% Store the class of the inputs
Iclass=class(Imoving);

% Convert the inputs to double
Imoving=double(Imoving);
Istatic=double(Istatic);
Imin=min(min(Istatic(:)),min(Imoving(:))); Imax=max(max(Istatic(:)),max(Istatic(:)));
Imoving=(Imoving-Imin)/(Imax-Imin);
Istatic=(Istatic-Imin)/(Imax-Imin);

function Ireg=Back2OldRange(Ireg,Iclass,Imin,Imax)
% Back to old image range
Ireg=Ireg*(Imax-Imin)+Imin;

% Set the class of output to input class
if(strcmpi(Iclass,'uint8')), Ireg=uint8(Ireg); end
if(strcmpi(Iclass,'uint16')), Ireg=uint16(Ireg); end
if(strcmpi(Iclass,'uint32')), Ireg=uint32(Ireg); end
if(strcmpi(Iclass,'int8')), Ireg=int8(Ireg); end
if(strcmpi(Iclass,'int16')), Ireg=int16(Ireg); end
if(strcmpi(Iclass,'int32')), Ireg=int32(Ireg); end
if(strcmpi(Iclass,'single')), Ireg=single(Ireg); end

function bsplinetranformimage()
data=getMyData();
data.RealGrid=data.Uniform+(data.Uniform-data.Grid);
data.Ireg=bspline_transform(data.RealGrid,data.Imoving,data.Spacing,3);
setMyData(data);

function showgrid()
data=getMyData();
if(ishandle(data.gridpointhandle))
    delete(data.gridpointhandle);
    delete(data.gridlinevhandle);
    delete(data.gridlinehhandle);
end
data.gridpointhandle=zeros(1,size(data.Grid,1)*size(data.Grid,2));
data.gridlinevhandle=zeros(1,size(data.Grid,1)*size(data.Grid,2));
data.gridlinehhandle=zeros(1,size(data.Grid,1)*size(data.Grid,2));

k=0;
hold(data.axesgrid_handle,'on');
IGrid=data.Grid+data.padding;
for i=1:size(IGrid,1),
    for j=1:size(IGrid,2),
        k=k+1;
        data.gridlinevhandle(k)=plot(data.axesgrid_handle,[IGrid(i,j,2) IGrid(i,min(j+1,size(IGrid,2)),2)],[ IGrid(i,j,1) IGrid(i,min(j+1,size(IGrid,2)),1)]);
        data.gridlinehhandle(k)=plot(data.axesgrid_handle,[IGrid(i,j,2) IGrid(min(i+1,size(IGrid,1)),j,2)],[IGrid(i,j,1) IGrid(min(i+1,size(IGrid,1)),j,1)]);
        data.gridpointhandle(k)=plot(data.axesgrid_handle,IGrid(i,j,2),IGrid(i,j,1),'go','MarkerSize',5);
        set(data.gridpointhandle(k),'ButtonDownFcn','manually_warp_images(''pointGridButtonDownFcn'',gcbo,[],guidata(gcbo))');

    end
end
hold(data.axesgrid_handle,'off');
setMyData(data);

function showslices()
data=getMyData();
x1=1+data.padding; x2=data.padding+size(data.Ireg,1);
y1=1+data.padding; y2=data.padding+size(data.Ireg,2);
if(size(data.Ireg,3)==1)
    data.Jreg=zeros([size(data.Ireg) 3]);
    data.Jreg(:,:,1)=data.Ireg; data.Jreg(:,:,2)=0.5*data.Ireg; 
else
    data.Jreg=data.Ireg;
end

switch(data.viewoption)
    case 0
        data.Itotal(x1:x2,y1:y2,:)=0.5*(data.Jreg+data.Jstatic);
    case 1
        diff=sum(abs(data.Ireg-data.Istatic),3); diff=diff./max(diff(:));
        data.Itotal(x1:x2,y1:y2,1)=diff;
        data.Itotal(x1:x2,y1:y2,2)=diff;
        data.Itotal(x1:x2,y1:y2,3)=diff;
    case 2
        diff=sum(abs(data.Ireg-data.Istatic),3); diff=diff*(255/max(diff(:)));
        data.Itotal(x1:x2,y1:y2,:)=ind2rgb(uint8(diff),jet(255));
    case 3
        data.Itotal(x1:x2,y1:y2,:)=data.Jreg;
    case 4
        data.Itotal(x1:x2,y1:y2,:)=data.Jstatic;
    case 5
        if(data.viewswap==0)
            data.Itotal(x1:x2,y1:y2,:)=data.Jreg;
        else
            data.Itotal(x1:x2,y1:y2,:)=data.Jstatic;
        end
    otherwise 
end
data.Itotal(data.Itotal<0)=0;
data.Itotal(data.Itotal>1)=1;
set(data.imshowhandle,'CData',data.Itotal); 
set(data.imshowpreview1handle,'CData',data.Ireg); 
drawnow

% --- Outputs from this function are returned to the command line.
function varargout = manually_warp_images_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function setMyData(data)
% Store data struct in figure
setappdata(gcf,'datareg',data);

function data=getMyData()
% Get data struct stored in figure
data=getappdata(gcf,'datareg');

function pointGridButtonDownFcn(hObject, eventdata, handles)
data=getMyData(); if(isempty(data)), return, end
data.pointselected=find(data.gridpointhandle==gcbo);
data.pointselectedhandle=gcbo;
set(data.pointselectedhandle, 'MarkerSize',8);
setMyData(data);

function TimerFcn(hObject, eventdata, handles)
data=getMyData();  if(isempty(data)), return, end
if(data.viewoption==5)
     data.viewswap=data.viewswap+1;
     if(data.viewswap>1), data.viewswap=0; end
     setMyData(data);
     showslices();
end

function cursor_position_update()
data=getMyData(); if(isempty(data)), return, end
p = get(data.axesgrid_handle, 'CurrentPoint');
data.mouse_position=[p(1, 1) p(1, 2)];
setMyData(data);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
if(ishandle(data.pointselectedhandle))
    set(data.pointselectedhandle, 'MarkerSize',5);
end
updateGrid();


function updateGrid()
data=getMyData();
[j,i]=ind2sub([size(data.Grid,2) size(data.Grid,1)],data.pointselected);
data.Grid(i,j,1)=data.mouse_position(2)-data.padding;
data.Grid(i,j,2)=data.mouse_position(1)-data.padding;
data.pointselectedhandle=[];
setMyData(data);
bsplinetranformimage();
showslices();
showgrid();

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cursor_position_update();
data=getMyData(); if(isempty(data)), return, end
if(~isempty(data.pointselectedhandle))
    set(data.pointselectedhandle, 'XData',data.mouse_position(1));
    set(data.pointselectedhandle, 'YData',data.mouse_position(2));
end


% --------------------------------------------------------------------
function menu_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_coloroverlap_Callback(hObject, eventdata, handles)
% hObject    handle to menu_coloroverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.viewoption=0;
setMyData(data);
set_menu_checks()
showslices();

% --------------------------------------------------------------------
function menu_diffgray_Callback(hObject, eventdata, handles)
% hObject    handle to menu_diffgray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.viewoption=1;
setMyData(data);
set_menu_checks()
showslices();


% --------------------------------------------------------------------
function menu_diffcolor_Callback(hObject, eventdata, handles)
% hObject    handle to menu_diffcolor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.viewoption=2;
setMyData(data);
set_menu_checks()
showslices();

% --------------------------------------------------------------------
function menu_movingimage_Callback(hObject, eventdata, handles)
% hObject    handle to menu_movingimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.viewoption=3;
setMyData(data);
set_menu_checks()
showslices();


% --------------------------------------------------------------------
function menu_staticimage_Callback(hObject, eventdata, handles)
% hObject    handle to menu_staticimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.viewoption=4;
setMyData(data);
set_menu_checks()
showslices();


% --------------------------------------------------------------------
function menu_swap_Callback(hObject, eventdata, handles)
% hObject    handle to menu_swap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData(); if(isempty(data)), return, end
data.viewoption=5;
data.viewswap=0;
setMyData(data);
set_menu_checks()
showslices();


% --------------------------------------------------------------------
function menu_registration_Callback(hObject, eventdata, handles)
% hObject    handle to menu_registration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_pixeldistance_Callback(hObject, eventdata, handles)
% hObject    handle to menu_pixeldistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
data.Similarity='sd';
setMyData(data);
set_menu_checks()

% --------------------------------------------------------------------
function menu_mutualinf_Callback(hObject, eventdata, handles)
% hObject    handle to menu_mutualinf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
data.Similarity='mi';
setMyData(data);
set_menu_checks()

% --------------------------------------------------------------------
function menu_startregistration_Callback(hObject, eventdata, handles)
% hObject    handle to menu_startregistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
Options.Spacing=data.Spacing;
Options.Grid=data.RealGrid;
Options.Similarity=data.Similarity;
Options.MaxRef=0;
Options.Registration='NonRigid';
[data.Ireg,data.RealGrid,data.Spacing] = image_registration(data.Imoving,data.Istatic,Options);
data.Grid=2*data.Uniform-data.RealGrid;
setMyData(data);
showslices();
showgrid();

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_grid_Callback(hObject, eventdata, handles)
% hObject    handle to menu_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_refinement_Callback(hObject, eventdata, handles)
% hObject    handle to menu_refinement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
[data.RealGrid,t]=refine_grid(data.RealGrid,data.Spacing,size(data.Imoving));
[data.Uniform,data.Spacing]=refine_grid(data.Uniform,data.Spacing,size(data.Imoving));
data.Grid=2*data.Uniform-data.RealGrid;
setMyData(data);
showslices();
showgrid();

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
data=getMyData();
% Stop the timer function
stop(data.GUItimer); delete(data.GUItimer); 
% Close the Window        
delete(hObject);


% --------------------------------------------------------------------
function menu_saveimage_Callback(hObject, eventdata, handles)
% hObject    handle to menu_saveimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile({'*.png';'*.jpg'}, 'Save Registerd Image as');
data=getMyData(); if(isempty(data)), return, end
imwrite(data.Ireg,[pathname filename]);

% --------------------------------------------------------------------
function menu_savematlab_Callback(hObject, eventdata, handles)
% hObject    handle to menu_savematlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uiputfile({'*.mat';'*.jpg'}, 'Save Matlab Variables as');
data=getMyData(); if(isempty(data)), return, end
Grid=data.RealGrid;
Spacing=data.Spacing;
Imoving=Back2OldRange(data.Imoving,data.Iclass,data.Imin,data.Imax);
Istatic=Back2OldRange(data.Istatic,data.Iclass,data.Imin,data.Imax);
[Ireg,Bx,By]=bspline_transform(Grid,Imoving,Spacing,3);
[Fx,Fy]=backwards2forwards(Bx,By);
save([filename pathname],'Grid','Spacing','Imoving','Istatic','Ireg','Bx','By','Fx','Fy');

% --------------------------------------------------------------------
function menu_resulttoworkspace_Callback(hObject, eventdata, handles)
% hObject    handle to menu_resulttoworkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
data=getMyData(); if(isempty(data)), return, end
Grid=data.RealGrid;
Spacing=data.Spacing;
Imoving=Back2OldRange(data.Imoving,data.Iclass,data.Imin,data.Imax);
Istatic=Back2OldRange(data.Istatic,data.Iclass,data.Imin,data.Imax);
[Ireg,Bx,By]=bspline_transform(Grid,Imoving,Spacing,3);

[Fx,Fy]=backwards2forwards(Bx,By);
assignin('base','Grid', Grid);
assignin('base','Spacing', Spacing);
assignin('base','Imoving', Imoving);
assignin('base','Istatic', Istatic);
assignin('base','Ireg',Ireg);
assignin('base','Bx', Bx);
assignin('base','By', By);
assignin('base','Fx', Fx);
assignin('base','Fy', Fy);

% --------------------------------------------------------------------
function menu_startaffine_Callback(hObject, eventdata, handles)
% hObject    handle to menu_startaffine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=getMyData();
Options.Spacing=data.Spacing;
Options.Grid=data.RealGrid;
Options.Similarity=data.Similarity;
Options.MaxRef=0;
Options.Registration='Affine';
[data.Ireg,data.RealGrid,data.Spacing] = image_registration(data.Imoving,data.Istatic,Options);
data.Grid=2*data.Uniform-data.RealGrid;
setMyData(data);
showslices();
showgrid();


function set_menu_checks()
data=getMyData();
set(data.handles.menu_coloroverlap,'Checked','off');
set(data.handles.menu_diffgray,'Checked','off');
set(data.handles.menu_diffcolor,'Checked','off');
set(data.handles.menu_movingimage,'Checked','off');
set(data.handles.menu_staticimage,'Checked','off');
set(data.handles.menu_swap,'Checked','off');
set(data.handles.menu_pixeldistance,'Checked','off');
set(data.handles.menu_mutualinf,'Checked','off');

switch(data.viewoption)
    case 0
        set(data.handles.menu_coloroverlap,'Checked','on');
    case 1
        set(data.handles.menu_diffgray,'Checked','on');
    case 2
        set(data.handles.menu_diffcolor,'Checked','on');
    case 3
        set(data.handles.menu_movingimage,'Checked','on');
    case 4
        set(data.handles.menu_staticimage,'Checked','on');
    case 5
        set(data.handles.menu_swap,'Checked','on');
    otherwise 
end

switch data.Similarity
    case 'mi'
        set(data.handles.menu_mutualinf,'Checked','on');
    case 'sd'
        set(data.handles.menu_pixeldistance,'Checked','on');
end
