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

% Last Modified by GUIDE v2.5 19-Apr-2016 14:30:05

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

Gui_Align_main(hObject, handles);





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



% change '/' to '\' due to the difference between python and matlab
failed_files_all = strrep(fileread('bad_files.txt'),'/','\');
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

cur_file = strrep(cur_file, '\Volumes\behavgenom_archive$\MaskedVideos', 'Z:\thecus');

avi_file = strrep(cur_file, '.hdf5', '.avi');

set(handles.edit1,'string',cur_file);
set(handles.text13,'string',avi_file);

system(['vlc "',avi_file,'"']);

%play

% f = figure('Position',[0 0 800 600]);
% a=actxcontrol('VideoLAN.VLCPlugin.2',[0,0,800,600]);
% a.playlist.add(avi_file);
% a.playlist.play();





