function varargout = Manual_Align(varargin)
% MANUAL_ALIGN MATLAB code for Manual_Align.fig
%      MANUAL_ALIGN, by itself, creates a new MANUAL_ALIGN or raises the existing
%      singleton*.
%
%      H = MANUAL_ALIGN returns the handle to a new MANUAL_ALIGN or the handle to
%      the existing singleton*.
%
%      MANUAL_ALIGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUAL_ALIGN.M with the given input arguments.
%
%      MANUAL_ALIGN('Property','Value',...) creates a new MANUAL_ALIGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Manual_Align_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Manual_Align_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Manual_Align

% Last Modified by GUIDE v2.5 27-Apr-2016 16:17:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Manual_Align_OpeningFcn, ...
    'gui_OutputFcn',  @Manual_Align_OutputFcn, ...
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


% --- Executes just before Manual_Align is made visible.
function Manual_Align_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Manual_Align (see VARARGIN)

% Choose default command line output for Manual_Align
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Manual_Align wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Manual_Align_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.checkbox1,'Value',0);
% set(handles.edit4,'Enable','off');
% set(handles.edit5,'Enable','off');
% set(handles.edit6,'Enable','off');

terminated =0;

set(handles.text18,'string','processing, please wait');
set(handles.pushbutton1,'Enable','on');
set(handles.pushbutton2,'Enable','on');
set(handles.pushbutton5,'Enable','on');
set(handles.checkbox1,'Enable','on');
set(handles.checkbox2,'Enable','on');
set(handles.checkbox3,'Enable','on');
set(handles.checkbox4,'Enable','on');
set(handles.checkbox5,'Enable','on');
drawnow()

% main function to do the alignment
[frame_diffs_d0,xyShift] = Gui_Align_main(hObject, handles, eventdata);

handles.frame_diffs_d0 = frame_diffs_d0;
handles.xyShift = xyShift;

guidata(hObject, handles);

set(handles.pushbutton6,'Enable','on');
set(handles.pushbutton7,'Enable','on');


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(gcf)




% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcf)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% clear figues
cla(handles.axes1)
cla(handles.axes2)
cla(handles.axes3)
cla(handles.text18)
%cla(handles.uitable1)
set(handles.uitable1,'data',[NaN,NaN;NaN,NaN]);

global terminated;
terminated =0;

if handles.edit2.Enable & handles.edit3.Enable
    
    read_file_list = get(handles.edit3,'String');
    
    try
        failed_files_all = strrep(fileread(read_file_list),'/','\');
    catch ME
        disp('The file list to read is not correct.');
        % change '/' to '\' due to the difference between python and matlab
        failed_files_all = strrep(fileread('bad_files.txt'),'/','\');
        % indicate otsuThr
        set(handles.edit3,'string','bad_files.txt');
    end
    % replace folder
    gap_sym = '\Volumes\behavgenom_archive$';
    
    ini_loc = strfind(failed_files_all,gap_sym);
    %ini_loc = regexp(failed_files_all,gap_sym);
    
    file_name = {};
    
    % restore file names to independent cell
    for ii = 1:numel(ini_loc)-1
        file_name = [file_name;failed_files_all(ini_loc(ii):ini_loc(ii+1)-2)];
    end
    file_name = [file_name;failed_files_all(ini_loc(numel(ini_loc)):end)];
    
    iif = str2num(get(handles.edit2,'string'));
    
    % set current file and result hdf5 file
    cur_file = strtrim(file_name{iif});
    
else
    cur_file = strtrim(get(handles.edit1,'String'));
    cur_file = strrep(cur_file,'/','\');
end

cur_file = strrep(cur_file, '\Volumes\behavgenom_archive$\', 'Z:\');

avi_file_temp = strrep(cur_file, '.hdf5', '.avi');
avi_file = strrep(avi_file_temp, 'MaskedVideos', 'thecus');

set(handles.edit1,'string',cur_file);
set(handles.text13,'string',avi_file);
try
    system(['vlc "',avi_file,'"']);
catch ME
    disp('The file cannot be read correctly.');
end

set(handles.pushbutton1,'Enable','on');
set(handles.text18,'string','watching the .avi video');

set(handles.pushbutton6,'Enable','off');
set(handles.pushbutton7,'Enable','off');

set(handles.checkbox1,'Enable','off');
set(handles.checkbox1,'Value',0);
set(handles.checkbox2,'Enable','off');
set(handles.checkbox2,'Value',0);
set(handles.checkbox3,'Enable','off');
set(handles.checkbox3,'Value',0);
set(handles.checkbox4,'Enable','off');
set(handles.checkbox4,'Value',0);
set(handles.checkbox5,'Enable','off');
set(handles.checkbox5,'Value',0);
set(handles.edit19,'Enable','off');
set(handles.edit20,'Enable','off');

set(handles.edit4,'Enable','off');
set(handles.edit5,'Enable','off');
set(handles.edit6,'Enable','off');
set(handles.edit7,'Enable','off');
set(handles.edit8,'Enable','off');
set(handles.edit9,'Enable','off');
set(handles.edit10,'Enable','off');
set(handles.edit11,'Enable','off');
set(handles.edit12,'Enable','off');
set(handles.edit13,'Enable','off');
set(handles.edit14,'Enable','off');
set(handles.edit15,'Enable','off');
set(handles.edit16,'Enable','off');
set(handles.edit17,'Enable','off');
set(handles.edit18,'Enable','off');






%play

% f = figure('Position',[0 0 800 600]);
% a=actxcontrol('VideoLAN.VLCPlugin.2',[0,0,800,600]);
% a.playlist.add(avi_file);
% a.playlist.play();


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global terminated;
terminated = 1;

set(handles.pushbutton1,'Enable','off');
set(handles.pushbutton2,'Enable','off');
set(handles.pushbutton3,'Enable','off');
set(handles.pushbutton5,'Enable','off');




function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

% handles.checkbx1 = get(hObject,'Value');

if handles.checkbox1.Value  %handles.checkbx1
    set(handles.edit4,'Enable','on');
    set(handles.edit5,'Enable','on');
    set(handles.edit6,'Enable','on');
    set(handles.edit20,'Enable','on');
    if strcmp(handles.pushbutton7.Enable,'on')
        set(handles.edit19,'Enable','on');
    end
else
    set(handles.edit4,'Enable','off');
    set(handles.edit5,'Enable','off');
    set(handles.edit6,'Enable','off');
    if ~(handles.checkbox1.Value | handles.checkbox2.Value |handles.checkbox3.Value...
            |handles.checkbox4.Value|handles.checkbox5.Value)
        set(handles.edit19,'Enable','off');
        set(handles.edit20,'Enable','off');
    end
end
% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
%handles.checkbx2 = get(hObject,'Value');

if handles.checkbox2.Value  %handles.checkbx2
    set(handles.edit7,'Enable','on');
    set(handles.edit8,'Enable','on');
    set(handles.edit9,'Enable','on');
    if strcmp(handles.pushbutton7.Enable,'on')
        set(handles.edit19,'Enable','on');
    end
    set(handles.edit20,'Enable','on');
else
    set(handles.edit7,'Enable','off');
    set(handles.edit8,'Enable','off');
    set(handles.edit9,'Enable','off');
    if ~(handles.checkbox1.Value | handles.checkbox2.Value |handles.checkbox3.Value...
            |handles.checkbox4.Value|handles.checkbox5.Value)
        set(handles.edit19,'Enable','off');
        set(handles.edit20,'Enable','off');
    end
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

if handles.checkbox3.Value  %handles.checkbx2
    set(handles.edit10,'Enable','on');
    set(handles.edit11,'Enable','on');
    set(handles.edit12,'Enable','on');
    if strcmp(handles.pushbutton7.Enable,'on')
        set(handles.edit19,'Enable','on');
    end
    set(handles.edit20,'Enable','on');
else
    set(handles.edit10,'Enable','off');
    set(handles.edit11,'Enable','off');
    set(handles.edit12,'Enable','off');
    if ~(handles.checkbox1.Value | handles.checkbox2.Value |handles.checkbox3.Value...
            |handles.checkbox4.Value|handles.checkbox5.Value)
        set(handles.edit19,'Enable','off');
        set(handles.edit20,'Enable','off');
    end
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
if handles.checkbox4.Value  %handles.checkbx2
    set(handles.edit13,'Enable','on');
    set(handles.edit14,'Enable','on');
    set(handles.edit15,'Enable','on');
    if strcmp(handles.pushbutton7.Enable,'on')
        set(handles.edit19,'Enable','on');
    end
    set(handles.edit20,'Enable','on');
else
    set(handles.edit13,'Enable','off');
    set(handles.edit14,'Enable','off');
    set(handles.edit15,'Enable','off');
    if ~(handles.checkbox1.Value | handles.checkbox2.Value |handles.checkbox3.Value...
            |handles.checkbox4.Value|handles.checkbox5.Value)
        set(handles.edit19,'Enable','off');
        set(handles.edit20,'Enable','off');
    end
end
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5

if handles.checkbox5.Value  %handles.checkbx2
    set(handles.edit16,'Enable','on');
    set(handles.edit17,'Enable','on');
    set(handles.edit18,'Enable','on');
    if strcmp(handles.pushbutton7.Enable,'on')
        set(handles.edit19,'Enable','on');
    end
    set(handles.edit20,'Enable','on');
else
    set(handles.edit16,'Enable','off');
    set(handles.edit17,'Enable','off');
    set(handles.edit18,'Enable','off');
    if ~(handles.checkbox1.Value | handles.checkbox2.Value |handles.checkbox3.Value...
            |handles.checkbox4.Value|handles.checkbox5.Value)
        set(handles.edit19,'Enable','off');
        set(handles.edit20,'Enable','off');
    end
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.checkbox1,'Value',0);
set(handles.edit4,'Enable','off');
set(handles.edit5,'Enable','off');
set(handles.edit6,'Enable','off');
set(handles.edit4,'Value',0);
set(handles.edit5,'Value',0);
set(handles.edit6,'Value',0);

set(handles.checkbox2,'Value',0);
set(handles.edit7,'Enable','off');
set(handles.edit8,'Enable','off');
set(handles.edit9,'Enable','off');
set(handles.edit7,'Value',0);
set(handles.edit8,'Value',0);
set(handles.edit9,'Value',0);

set(handles.checkbox3,'Value',0);
set(handles.edit10,'Enable','off');
set(handles.edit11,'Enable','off');
set(handles.edit12,'Enable','off');
set(handles.edit10,'Value',0);
set(handles.edit11,'Value',0);
set(handles.edit12,'Value',0);

set(handles.checkbox4,'Value',0);
set(handles.edit13,'Enable','off');
set(handles.edit14,'Enable','off');
set(handles.edit15,'Enable','off');
set(handles.edit13,'Value',0);
set(handles.edit14,'Value',0);
set(handles.edit15,'Value',0);

set(handles.checkbox5,'Value',0);
set(handles.edit16,'Enable','off');
set(handles.edit17,'Enable','off');
set(handles.edit18,'Enable','off');
set(handles.edit16,'Value',0);
set(handles.edit17,'Value',0);
set(handles.edit18,'Value',0);

set(handles.pushbutton6,'Enable','off');
set(handles.pushbutton7,'Enable','off');
set(handles.edit20,'Enable','off');


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Gui_Align_rerun(hObject, handles);



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
if handles.checkbox7.Value  %handles.checkbx2
    set(handles.edit1,'Enable','on');
    set(handles.edit2,'Enable','off');
    set(handles.edit3,'Enable','off');
else
    set(handles.edit1,'Enable','off');
    set(handles.edit2,'Enable','on');
    set(handles.edit3,'Enable','on');
end
% Update handles structure
guidata(hObject, handles);




function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if strcmp( get(handles.figure1,'selectionType') , 'normal')
%     disp('Left Click')
% end
if strcmp( get(handles.figure1,'selectionType') , 'open')
    disp('Left Double Click')
	%set(handles.text27,'string',['(',num2str(eventdata.IntersectionPoint(1)), num2str(eventdata.IntersectionPoint(2)),')']);
    set(handles.text27,'string',['(',sprintf('%5.2f',eventdata.IntersectionPoint(1)),',',sprintf('%5.2f',eventdata.IntersectionPoint(2)),')']);
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
set (gcf, 'WindowButtonMotionFcn', @mouseMove);


function mouseMove (hObject, eventdata, handles)
C = get (gca, 'CurrentPoint');
x_cord = num2str(C(1,1));
y_cord = num2str(C(1,2));
title(gca, [ '\fontsize{9}(X,Y) = (',x_cord, ', ',y_cord, ')'],'Position', [0.4 1]);



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit4.
function edit4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Left Click Edit')
cord_xy = get(handles.text27,'string');
[cord_x, cord_y] = strtok(cord_xy, '(,');

%set(handles.text27,'string',['(',num2str(eventdata.IntersectionPoint(1)), num2str(eventdata.IntersectionPoint(2)),')']);
set(handles.text27,'string',['(',sprintf('%5.2f',eventdata.IntersectionPoint(1)),',',sprintf('%5.2f',eventdata.IntersectionPoint(2)),')']);



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
pop_cont = get(hObject,'Value');

disp('Left Click Edit')
cord_xy = get(handles.text27,'string');
[cord_x, rem] = strtok(cord_xy, '(,');
[cord_y, rem2] = strtok(rem, ',)');


switch pop_cont
    case 2
        set(handles.edit4,'string',cord_x);
    case 3
        set(handles.edit5,'string',cord_x);
        set(handles.edit6,'string',cord_y);
    case 4
        set(handles.edit7,'string',cord_x);
    case 5
        set(handles.edit8,'string',cord_x);
        set(handles.edit9,'string',cord_y);
    case 6
        set(handles.edit10,'string',cord_x);
    case 7
        set(handles.edit11,'string',cord_x);
        set(handles.edit12,'string',cord_y);
    case 8
        set(handles.edit13,'string',cord_x);
    case 9
        set(handles.edit14,'string',cord_x);
        set(handles.edit15,'string',cord_y);
    case 10
        set(handles.edit16,'string',cord_x);
    case 11
        set(handles.edit17,'string',cord_x);
        set(handles.edit18,'string',cord_y);
    otherwise
        set(handles.text18,'string','nothing is copied to text box');
end



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
