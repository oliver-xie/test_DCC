function varargout = simtb_GUI_matrix_gen_B(varargin)
% SIMTB_GUI_MATRIX_GEN_B M-file for simtb_GUI_matrix_gen_B.fig
%      SIMTB_GUI_MATRIX_GEN_B, by itself, creates a new SIMTB_GUI_MATRIX_GEN_B or raises the existing
%      singleton*.
%
%      H = SIMTB_GUI_MATRIX_GEN_B returns the handle to a new SIMTB_GUI_MATRIX_GEN_B or the handle to
%      the existing singleton*.
%
%      SIMTB_GUI_MATRIX_GEN_B('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMTB_GUI_MATRIX_GEN_B.M with the given input arguments.
%
%      SIMTB_GUI_MATRIX_GEN_B('Property','Value',...) creates a new SIMTB_GUI_MATRIX_GEN_B or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simtb_GUI_matrix_gen_B_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to simtb_GUI_matrix_gen_B_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help simtb_GUI_matrix_gen_B

% Last Modified by GUIDE v2.5 15-Dec-2010 15:41:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @simtb_GUI_matrix_gen_B_OpeningFcn, ...
    'gui_OutputFcn',  @simtb_GUI_matrix_gen_B_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);

% Note: the first argument is not a callback function
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before simtb_GUI_matrix_gen_B is made visible.
function simtb_GUI_matrix_gen_B_OpeningFcn(hObject, eventdata, handles, varargin)
set(handles.figure_gen, 'Color', [0 0 0])
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simtb_GUI_matrix_gen_B (see VARARGIN)
handles.fgColor = [0.6784    0.9216    1.0000];
handles.bgColor = [.5 .5 .5];

if nargin < 2
    handles.sP.M=10;
    handles.sP.nC = 10;
else
    axxlabel = varargin{2};
    axylabel = varargin{3};
    cnames = varargin{4};
    rnames = varargin{5};
    instruction = varargin{6};
    default = varargin{7};
    name = varargin{8};
  %  paneltitle = varargin{9};
    paramname = varargin{10};
   
    
    columnformat = cell(1,length(cnames));
    [columnformat{:}] = deal('short g');
    set(handles.uitable1, 'ColumnFormat', columnformat);
   
    handles.paramname = paramname;
   % set(handles.panel_title, 'String', paneltitle);
    set(handles.figure_gen,'Name', name);
    set(handles.instruction, 'Max',2, 'ForeGroundcolor', handles.fgColor, 'Fontsize', 10, 'String', instruction, 'HorizontalAlignment', 'left');
    set(handles.xlabel,'String', axxlabel);
    handles.ylabel =  axylabel;
    set(handles.uitable1,'ColumnName',cnames,'RowName',rnames);
    set(handles.uitable1, 'Data', default);
    set(handles.uitable1,'Visible','On');
    [handles.row, handles.col] = size(default);
    handles.default = default;
end

handles.data = ones(handles.row, handles.col);


%% adjust the position of the table based on the extent
%--------------------------------------------------------------------------
set(handles.figure_gen, 'Visible', 'on')
TPOS = get(handles.uitable1, 'Position');
LEFT = TPOS(1);
BOTTOM = TPOS(2);
WIDTH = TPOS(3);
HEIGHT = TPOS(4);
TOP = BOTTOM+HEIGHT;

TEXTENT = get(handles.uitable1, 'Extent');
fixedW = TEXTENT(3);
fixedH = TEXTENT(4);

minH = 3; % space for row/column headers
minW = 6;

if WIDTH > (fixedW+minW)
    WIDTH =  fixedW + minW;
end
if HEIGHT > (fixedH+minH)
    HEIGHT = fixedH + minH;
end
BOTTOM = TOP-HEIGHT;

set(handles.uitable1, 'Position', [LEFT ,BOTTOM, WIDTH, HEIGHT])
set(handles.figure_gen, 'Units', 'Normalized')
set(handles.uitable1, 'Units', 'Normalized')
TPOS = get(handles.uitable1, 'Position');
LEFT = TPOS(1);
BOTTOM = TPOS(2);
WIDTH = TPOS(3);
HEIGHT = TPOS(4);
set(handles.uitable1, 'Position', [0.5-WIDTH/2 ,BOTTOM, WIDTH, HEIGHT])
axPos = get(handles.axesY, 'Position');
set(handles.axesY, 'Position', [0.46-WIDTH/2, BOTTOM-0, axPos(3), HEIGHT], 'XLim', [-1 1], 'YLim', [-1 1])
T = text(0, 0, handles.ylabel, 'Rotation', 90, 'Fontsize', 8,'FontWeight', 'bold', ... 
    'HandleVisibility', 'on','Color', 'white', 'HorizontalAlignment', 'center');
uistack(T, 'top')
axis off
%-------------------------------------------------------------------------
% Update handles structure
guidata(hObject, handles);
uiwait;

% --- Outputs from this function are returned to the command line.
function varargout = simtb_GUI_matrix_gen_B_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.data = get(handles.uitable1, 'data');
varargout{1} = handles.data;
delete(handles.figure_gen);


function instruction_Callback(hObject, eventdata, handles)
% hObject    handle to instruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of instruction as text
%        str2double(get(hObject,'String')) returns contents of instruction as a double


% --- Executes during object creation, after setting all properties.
function instruction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to instruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu.
function popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu
% Determine the selected data set.
str = get(hObject, 'String');
val = get(hObject,'Value');
% Set current data to the selected data set.
% Close other buttons
switch str{val};
    case 'All Ones' % Same value for all elements.
        set(handles.tag,'Visible','Off');
        set(handles.value,'Visible','Off');
        temp = ones(handles.row, handles.col);
        set(handles.uitable1, 'data', temp);
        set(handles.uitable1,'Visible','On');

    case 'All Zeros' %Matrix elements are Uniform Distributed.
        set(handles.tag,'Visible','Off');
        set(handles.value,'Visible','Off');
        temp = zeros(handles.row, handles.col);
        set(handles.uitable1, 'data', temp);
        set(handles.uitable1,'Visible','On');

    case 'Random Zero/Ones' % Matrix elements are Randomly Distributed.
        set(handles.value,'Visible','On');
        set(handles.tag,'Visible','On');
end

% Save the handles structure.
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function value_Callback(hObject, eventdata, handles)
% hObject    handle to value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value as text
%        str2double(get(hObject,'String')) returns contents of value as a double
propotion = str2num(get(hObject,'String'));
temp = zeros(handles.row, handles.col);
temp(find(rand(handles.row, handles.col) < propotion )) = 1;
set(handles.uitable1, 'data', temp);
set(handles.uitable1,'Visible','On');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xlabel_Callback(hObject, eventdata, handles)
% hObject    handle to xlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlabel as text
%        str2double(get(hObject,'String')) returns contents of xlabel as a double


% --- Executes during object creation, after setting all properties.
function xlabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.popupmenu,'Value',1);
set(handles.tag,'Visible','Off');
set(handles.value,'Visible','Off');
set(handles.value,'Visible','Off');
set(handles.value,'String','');
set(handles.uitable1, 'data', handles.default);
set(handles.uitable1,'Visible','On');
set(hObject,'BackgroundColor',handles.fgColor);
set(handles.Done,'BackgroundColor', handles.bgColor);
guidata(hObject, handles);

% --- Executes on button press in Done.
function Done_Callback(hObject, eventdata, handles)
% hObject    handle to Done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject,'BackgroundColor',handles.fgColor);
set(handles.reset,'BackgroundColor', handles.bgColor);
uiresume;
guidata(hObject, handles);


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Title = 'Instruction';
% Icon = 'none';
% a = 'All Ones: All the element values are ones. ';
% b = 'All Zeros: All the element values are zeros.'; 
% c = 'Random Zero/Ones:  Element values are randomly ones/zeros with ones show at the proportion of the filling value.';
% Message = strvcat(a,b,c);
% h = msgbox(Message,Title,Icon);
simtb_GUI_format_paramhelp(handles.paramname);
