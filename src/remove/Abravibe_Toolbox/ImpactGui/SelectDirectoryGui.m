function varargout = SelectDirectoryGui(varargin)
% SELECTDIRECTORYGUI MATLAB code for SelectDirectoryGui.fig
%      SELECTDIRECTORYGUI, by itself, creates a new SELECTDIRECTORYGUI or raises the existing
%      singleton*.
%
%      H = SELECTDIRECTORYGUI returns the handle to a new SELECTDIRECTORYGUI or the handle to
%      the existing singleton*.
%
%      SELECTDIRECTORYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTDIRECTORYGUI.M with the given input arguments.
%
%      SELECTDIRECTORYGUI('Property','Value',...) creates a new SELECTDIRECTORYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectDirectoryGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectDirectoryGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectDirectoryGui

% Last Modified by GUIDE v2.5 13-Dec-2015 17:39:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectDirectoryGui_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectDirectoryGui_OutputFcn, ...
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


% --- Executes just before SelectDirectoryGui is made visible.
function SelectDirectoryGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectDirectoryGui (see VARARGIN)

% Choose default command line output for SelectDirectoryGui
handles.output = hObject;

% Fill GUI
handles.InputDir=pwd;
set(handles.InputDirText,'String',handles.InputDir);
handles=UpdateTable(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectDirectoryGui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SelectDirectoryGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles,'DirectoryList')
    % Build full paths for output of directories, as to this point only the
    % directory names are stored.
    for n=1:length(handles.SelectedDirectories)
        handles.DirectoryList{handles.SelectedDirectories(n)}=fullfile(handles.InputDir,handles.DirectoryList{handles.SelectedDirectories(n)});
    end
    handles.output=handles.DirectoryList(handles.SelectedDirectories);
    varargout{1} = handles.output;
    delete(handles.figure1)
else
    varargout{1} = [];
end


% --- Executes on button press in InputDirButton.
function InputDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to InputDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SelectedDir=uigetdir();
if SelectedDir ~= 0
    handles.InputDir=SelectedDir;
    set(handles.InputDirText,'String',handles.InputDir);
    handles.InputDir=SelectedDir;
    handles=UpdateTable(handles);
end
guidata(hObject,handles)

% --- Executes on button press in SelectButton.
function SelectButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;

% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;
delete(handles.figure1);

% --- Executes when selected cell(s) is changed in DirectoryTable.
function DirectoryTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to DirectoryTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
selection = eventdata.Indices(:,1);
handles.SelectedDirectories=selection;
guidata(hObject, handles);


%==========================================================================
% User defined functions
function  handles = UpdateTable(handles)
D=dir(fullfile(handles.InputDir,'*.'));
DirectoryNames={D.name};
DirectoryNames={DirectoryNames{3:end}}';     % Delete . and ..
set(handles.DirectoryTable,'Data',DirectoryNames)
handles.DirectoryList=DirectoryNames;
handles.SelectedDirectories=[1:length(DirectoryNames)]';
