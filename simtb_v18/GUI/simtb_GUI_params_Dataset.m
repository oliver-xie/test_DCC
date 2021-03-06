function varargout = simtb_GUI_params_Dataset(varargin)
% SIMTB_GUI_PARAMS_DATASET M-file for simtb_GUI_params_Dataset.fig
%      SIMTB_GUI_PARAMS_DATASET, by itself, creates SM_present new SIMTB_GUI_PARAMS_DATASET or raises the existing
%      singleton*.
%
%      H = SIMTB_GUI_PARAMS_DATASET returns the handle to SM_present new SIMTB_GUI_PARAMS_DATASET or the handle to
%      the existing singleton*.
%
%      SIMTB_GUI_PARAMS_DATASET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIMTB_GUI_PARAMS_DATASET.M with the given input arguments.
%
%      SIMTB_GUI_PARAMS_DATASET('Property','Value',...) creates SM_present new SIMTB_GUI_PARAMS_DATASET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before simtb_GUI_params_Dataset_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to simtb_GUI_params_Dataset_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help simtb_GUI_params_Dataset

% Last Modified by GUIDE v2.5 27-Nov-2010 20:37:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simtb_GUI_params_Dataset_OpeningFcn, ...
                   'gui_OutputFcn',  @simtb_GUI_params_Dataset_OutputFcn, ...
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


% --- Executes just before simtb_GUI_params_Dataset is made visible.
function simtb_GUI_params_Dataset_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in SM_present future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to simtb_GUI_params_Dataset (see VARARGIN)
handles.fgColor = [0.6784    0.9216    1.0000];
handles.bgColor = [.5 .5 .5];

set(handles.figure1, 'Name', 'Parameter Selection: Step 5');
set(handles.text1,  'HorizontalAlignment', 'left');
set(handles.text2,  'HorizontalAlignment', 'left');
set(handles.text3,  'HorizontalAlignment', 'left');
set(handles.text4,  'HorizontalAlignment', 'left');
set(handles.text5,  'HorizontalAlignment', 'left');
set(handles.text6,  'HorizontalAlignment', 'left');
set(handles.text7,  'HorizontalAlignment', 'left');
set(handles.text8,  'HorizontalAlignment', 'left');
set(handles.text9,  'HorizontalAlignment', 'left');

sP = varargin{1};
handles.sP=sP;

set(handles.D_motion_TRANSmax, 'String', num2str(handles.sP.D_motion_TRANSmax))
set(handles.D_motion_ROTmax, 'String', num2str(handles.sP.D_motion_ROTmax))

if handles.sP.D_TT_FLAG
    set(handles.D_TT_FLAG_yes,'BackgroundColor',handles.fgColor);
    set(handles.D_TT_FLAG_no,'BackgroundColor',handles.bgColor);
else
    set(handles.D_TT_FLAG_no,'BackgroundColor',handles.fgColor);
    set(handles.D_TT_FLAG_yes,'BackgroundColor',handles.bgColor);
end


if handles.sP.D_noise_FLAG
    set(handles.D_noise_FLAG_yes,'BackgroundColor',handles.fgColor);
    set(handles.D_noise_FLAG_no,'BackgroundColor',handles.bgColor);
else
    set(handles.D_noise_FLAG_no,'BackgroundColor',handles.fgColor);
    set(handles.D_noise_FLAG_yes,'BackgroundColor',handles.bgColor);
end

if handles.sP.D_motion_FLAG
    set(handles.D_motion_FLAG_yes,'BackgroundColor',handles.fgColor);
    set(handles.D_motion_FLAG_no,'BackgroundColor',handles.bgColor);
else
    set(handles.D_motion_FLAG_no,'BackgroundColor',handles.fgColor);
    set(handles.D_motion_FLAG_yes,'BackgroundColor',handles.bgColor);
end

set(handles.text1, 'String', simtb_GUI_format_paramlabel('D_baseline'), 'HorizontalAlignment', 'left');
set(handles.text2, 'String', simtb_GUI_format_paramlabel('D_TT_FLAG'), 'HorizontalAlignment', 'left');
set(handles.text3, 'String', simtb_GUI_format_paramlabel('D_pSC'), 'HorizontalAlignment', 'left');
set(handles.text4, 'String', simtb_GUI_format_paramlabel('D_noise_FLAG'), 'HorizontalAlignment', 'left');
set(handles.text5, 'String', simtb_GUI_format_paramlabel('D_CNR'), 'HorizontalAlignment', 'left');
set(handles.text6, 'String', simtb_GUI_format_paramlabel('D_motion_FLAG'), 'HorizontalAlignment', 'left');
set(handles.text7, 'String', simtb_GUI_format_paramlabel('D_motion_TRANSmax'), 'HorizontalAlignment', 'left');
set(handles.text8, 'String', simtb_GUI_format_paramlabel('D_motion_ROTmax'), 'HorizontalAlignment', 'left');
set(handles.text9, 'String', simtb_GUI_format_paramlabel('D_motion_deviates'), 'HorizontalAlignment', 'left');


guidata(hObject, handles);
handles.output = handles.sP;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes simtb_GUI_params_Dataset wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = simtb_GUI_params_Dataset_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in SM_present future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in D_baseline.
function D_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to D_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlabel = 'Subject Index';
ylabel = '';  %ylabel = 'Dataset';
cnames={};
rnames={};

[m,n] = size(handles.sP.D_baseline);
default = handles.sP.D_baseline;  
mxsize.row=1;
mxsize.col=handles.sP.M;

if m~=mxsize.row || n~=mxsize.col
   m=mxsize.row;
   n=mxsize.col;
   default = 800*ones(m,n);
end

for i= 1:m
  rnames{i} = '';
end

for i= 1:n
%   cnames{i} = num2str(i);
  cnames{i} = sprintf('      %3d',i);  % 9 character column
end

instruction = 'Fill in a baseline signal intensity for each subject. To get started, use the drop-down menu at left and select a method to initialize values.  Then, update the values by hand.  When you are satisfied, hit "Done". ';
panelname = 'Parameter Selection: Step 5a';
paneltitle = 'Baseline Intensity';
handles.sP.D_baseline = simtb_GUI_matrix_gen(2, xlabel, ylabel, cnames, rnames,instruction, default,panelname, paneltitle, 'D_baseline');
guidata(hObject, handles); 
set(hObject,'BackgroundColor',handles.fgColor);


% --- Executes on button press in D_TT_FLAG_yes.
function D_TT_FLAG_yes_Callback(hObject, eventdata, handles)
% hObject    handle to D_TT_FLAG_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sP.D_TT_FLAG = 1;
[TT_count, TT_levels] = simtb_countTT; 
handles.sP.D_TT_levels = TT_levels;

guidata(hObject, handles);
set(hObject,'BackgroundColor',handles.fgColor);
set(handles.D_TT_FLAG_no,'BackgroundColor',handles.bgColor);

Title = 'Using Defaults';
Icon = 'none';
[TT_count, TT_levels] = simtb_countTT;
TTstring = sprintf('%0.1f, ', TT_levels); TTstring(end) = ''; TTstring(end) = '';
Message = ['Using default Tissue Type levels [' TTstring '].'... 
' To customize these levels update the parameter structure field ''D_TT_levels'' or edit the function ''simtb_countTT''.'];
h = msgbox(Message,Title,Icon);


% --- Executes on button press in D_TT_FLAG_no.
function D_TT_FLAG_no_Callback(hObject, eventdata, handles)
% hObject    handle to D_TT_FLAG_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sP.D_TT_FLAG = 0;
guidata(hObject, handles);
set(hObject,'BackgroundColor',handles.fgColor);
set(handles.D_TT_FLAG_yes,'BackgroundColor',handles.bgColor);


% --- Executes on button press in D_pSC.
function D_pSC_Callback(hObject, eventdata, handles)
% hObject    handle to D_pSC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlabel = 'Source ID (Component Index)';
ylabel = 'Subject Index';
cnames={};
rnames={};

[m,n] = size(handles.sP.D_pSC);
default = handles.sP.D_pSC;  
mxsize.row=handles.sP.M;
mxsize.col=handles.sP.nC;

if m~=mxsize.row || n~=mxsize.col
   m=mxsize.row;
   n=mxsize.col;
   default = 3*ones(m,n);
end

for i= 1:m
%   rnames{i} = num2str(i);
  rnames{i} = sprintf('      %3d',i);  % 9 character column
end

for i= 1:n
%   cnames{i} = strcat(num2str(handles.sP.SM_source_ID(i)), ' (', num2str(i), ')');
  %   cnames{i} = num2str(i);
  cnames{i} = sprintf('%3d (%3d)',handles.sP.SM_source_ID(i),i);  % 9 character column
end

instruction = 'Fill in the  matrix of percent signal changes for each subject and component. To get started, use the drop-down menu at left and select a method to initialize values.  Then, update the values by hand.  When you are satisfied, hit "Done". ';
panelname = 'Parameter Selection: Step 5b';
paneltitle = 'Percent Signal Change from Baseline';
handles.sP.D_pSC = simtb_GUI_matrix_gen(2, xlabel, ylabel, cnames, rnames,instruction, default,panelname, paneltitle, 'D_pSC');
guidata(hObject, handles); 
set(hObject,'BackgroundColor',handles.fgColor);


% --- Executes on button press in D_noise_FLAG_yes.
function D_noise_FLAG_yes_Callback(hObject, eventdata, handles)
% hObject    handle to D_noise_FLAG_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sP.D_noise_FLAG  = 1;
guidata(hObject, handles);
set(hObject,'BackgroundColor',handles.fgColor);
set(handles.D_noise_FLAG_no,'BackgroundColor',handles.bgColor);

% --- Executes on button press in D_noise_FLAG_no.
function D_noise_FLAG_no_Callback(hObject, eventdata, handles)
% hObject    handle to D_noise_FLAG_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sP.D_noise_FLAG  = 0;
guidata(hObject, handles);
set(hObject,'BackgroundColor',handles.fgColor);
set(handles.D_noise_FLAG_yes,'BackgroundColor',handles.bgColor);


% --- Executes on button press in D_CNR.
function D_CNR_Callback(hObject, eventdata, handles)
% hObject    handle to D_CNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlabel = 'Subject Index';
ylabel = '';
cnames={};
rnames={};

[m,n] = size(handles.sP.D_CNR);
default = handles.sP.D_CNR;  
mxsize.row=1;
mxsize.col=handles.sP.M;

if m~=mxsize.row || n~=mxsize.col
   m=mxsize.row;
   n=mxsize.col;
   default = zeros(m,n);
end

for i= 1:m
  rnames{i} = '';
end

for i= 1:n
%   cnames{i} =  num2str(i);
  cnames{i} = sprintf('      %3d',i);  % 9 character column
end

instruction = 'Fill in the vector of contrast-to-noise ratio for each subject. To get started, use the drop-down menu at left and select a method to initialize values.  Then, update the values by hand.  When you are satisfied, hit "Done". ';
panelname = 'Parameter Selection: Step 5c';
paneltitle = 'Contrast-to-Noise Ratio';
handles.sP.D_CNR = simtb_GUI_matrix_gen(2,xlabel, ylabel, cnames, rnames,instruction, default,panelname, paneltitle, 'D_CNR');
guidata(hObject, handles); 
set(hObject,'BackgroundColor',handles.fgColor);


% --- Executes on button press in D_motion_FLAG_yes.
function D_motion_FLAG_yes_Callback(hObject, eventdata, handles)
% hObject    handle to D_motion_FLAG_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sP.D_motion_FLAG  = 1;
guidata(hObject, handles);
set(hObject,'BackgroundColor',handles.fgColor);
set(handles.D_motion_FLAG_no,'BackgroundColor',handles.bgColor);


% --- Executes on button press in D_motion_FLAG_no.
function D_motion_FLAG_no_Callback(hObject, eventdata, handles)
% hObject    handle to D_motion_FLAG_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sP.D_motion_FLAG  = 0;
guidata(hObject, handles);
set(hObject,'BackgroundColor',handles.fgColor);
set(handles.D_motion_FLAG_yes,'BackgroundColor',handles.bgColor);


function D_motion_TRANSmax_Callback(hObject, eventdata, handles)
% hObject    handle to D_motion_TRANSmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D_motion_TRANSmax as text
%        str2double(get(hObject,'String')) returns contents of D_motion_TRANSmax as a double
handles.sP.D_motion_TRANSmax = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function D_motion_TRANSmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D_motion_TRANSmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function D_motion_ROTmax_Callback(hObject, eventdata, handles)
% hObject    handle to D_motion_ROTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D_motion_ROTmax as text
%        str2double(get(hObject,'String')) returns contents of D_motion_ROTmax as a double
handles.sP.D_motion_ROTmax = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function D_motion_ROTmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D_motion_ROTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in D_motion_deviates.
function D_motion_deviates_Callback(hObject, eventdata, handles)
% hObject    handle to D_motion_deviates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
xlabel = ' ';
ylabel = 'Subject Index';
cnames={};
rnames={};

[m,n] = size(handles.sP.D_motion_deviates);
default = handles.sP.D_motion_deviates;  
mxsize.row=handles.sP.M;
mxsize.col=3;

if m~=mxsize.row || n~=mxsize.col
   m=mxsize.row;
   n=mxsize.col;
   default = zeros(m,n);
end

for i= 1:m
%   rnames{i} = num2str(i);
  rnames{i} = sprintf('      %3d',i);  % 9 character column
end

%cnames{1} = sprintf('x-translation');
%cnames{2} = sprintf('y-translation');
%cnames{3} = sprintf('rotation');
cnames{1} = '  x-trans';
cnames{2} = '  y-trans';
cnames{3} = ' rotation'; % 9 character column

instruction = 'Fill in the matrix of motion deviates, in fractional units relative to the maximum motion for each subject. To get started, use the drop-down menu at left and select a method to initialize values.  Then, update the values by hand.  When you are satisfied, hit "Done". ';
panelname = 'Parameter Selection: Step 5d';
paneltitle = 'Motion Deviates';
handles.sP.D_motion_deviates = simtb_GUI_matrix_gen(2, xlabel, ylabel, cnames, rnames,instruction, default, panelname, paneltitle, 'D_motion_deviates');

%% Value for component presence must be either 0 or 1
if any(handles.sP.D_motion_deviates(:) < 0) || any(handles.sP.D_motion_deviates(:) > 1)
   Title ='';
   Icon ='error';
   Message = 'ERROR: Motion deviates must be between 0 and 1, representing a fraction of max motion. Re-enter ''D_motion_deviates'' value.';
   msgbox(Message,Title,Icon)
   set(hObject,'BackgroundColor','red');
else
   set(hObject,'BackgroundColor',handles.fgColor);
end
%%

guidata(hObject, handles); 

% --- Executes on button press in Back.
function Back_Callback(hObject, eventdata, handles)
% hObject    handle to Back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf);
simtb_GUI_params_TC_Basic(handles.sP);


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in SM_present future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[errorflag, Message] = simtb_checkparams(handles.sP, '6');
if errorflag
   Title = 'Error';
   Icon = 'error';
   h = msgbox(Message,Title,Icon);
else
         
  close(gcf);
  handles.sP = simtb_GUI_params_Final(handles.sP);             
end


% --- Executes on button press in help_D_baseline.
function help_D_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_baseline');

% --- Executes on button press in help_D_TT_FLAG.
function help_D_TT_FLAG_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_TT_FLAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_TT_FLAG');

% --- Executes on button press in help_D_pSC.
function help_D_pSC_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_pSC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_pSC');

% --- Executes on button press in help_D_noise_FLAG.
function help_D_noise_FLAG_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_noise_FLAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_noise_FLAG');

% --- Executes on button press in help_D_CNR.
function help_D_CNR_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_CNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_CNR');

% --- Executes on button press in help_D_motion_FLAG.
function help_D_motion_FLAG_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_motion_FLAG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_motion_FLAG');

% --- Executes on button press in help_D_motion_TRANSmax.
function help_D_motion_TRANSmax_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_motion_TRANSmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_motion_TRANSmax');

% --- Executes on button press in help_D_motion_ROTmax.
function help_D_motion_ROTmax_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_motion_ROTmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_motion_ROTmax');

% --- Executes on button press in help_D_motion_deviates.
function help_D_motion_deviates_Callback(hObject, eventdata, handles)
% hObject    handle to help_D_motion_deviates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
simtb_GUI_format_paramhelp('D_motion_deviates');
