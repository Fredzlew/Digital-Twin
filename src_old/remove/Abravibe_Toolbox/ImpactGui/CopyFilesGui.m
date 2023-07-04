function varargout = CopyFilesGui(varargin)
% COPYFILESGUI MATLAB code for CopyFilesGui.fig
%      COPYFILESGUI, by itself, creates a new COPYFILESGUI or raises the existing
%      singleton*.
%
%      H = COPYFILESGUI returns the handle to a new COPYFILESGUI or the handle to
%      the existing singleton*.
%
%      COPYFILESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COPYFILESGUI.M with the given input arguments.
%
%      COPYFILESGUI('Property','Value',...) creates a new COPYFILESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CopyFilesGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CopyFilesGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2009-2014 by Anders Brandt
% Email: abra@iti.sdu.dk
% Version: 1.0 2014-07-15
% This file is part of ABRAVIBE Toolbox for NVA


% Last Modified by GUIDE v2.5 18-Dec-2015 11:50:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CopyFilesGui_OpeningFcn, ...
                   'gui_OutputFcn',  @CopyFilesGui_OutputFcn, ...
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


% --- Executes just before CopyFilesGui is made visible.
function CopyFilesGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CopyFilesGui (see VARARGIN)

% Choose default command line output for CopyFilesGui
handles.output = '';

n=nargin;
if nargin >= 3
    a = varargin{2};
    b = varargin{3};
else
    errordlg('Wrong number of inputs in MergeFilesGui! Must be >= 3!');
end
handles.InputDirectory=a;
handles.OutputDirectory=b;
UpdateFileTable(hObject,handles);
set(handles.DirectoryText,'String',handles.OutputDirectory)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CopyFilesGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CopyFilesGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SelectFilesButton.
function SelectFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Files,Path]=uigetfile('*.imptime','Select Files to merge into output directory!',...
    handles.OutputDirectory,'MultiSelect','on');
% Check that files have not been selected that are already in the current
% output directory
if strcmp(Path,strcat(handles.OutputDirectory,filesep))
    errordlg('You cannot select files that are in the output directory!')
else
    if ~iscell(Files)
        if Files ~= 0
            Files={Files};
        end
    end
    if iscell(Files)
        for n = 1:length(Files)
            MergeFile(handles.OutputDirectory,fullfile(Path,Files{n}));
        end
    end
    UpdateFileTable(hObject,handles);
    guidata(hObject,handles)
end

% --- Executes on button press in DeleteFilesButton.
function DeleteFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check that one or more files have been selected in the table
if isfield(handles,'selection') & ~isempty(handles.selection);
    Files={handles.FileTable.Data{handles.selection}}
    if iscell(Files)
        s=questdlg('Warning! Files will be permanently deleted! Are you sure?');
        if strcmp(upper(s),'YES')
            for n = 1:length(Files)
                FileName=fullfile(handles.OutputDirectory,Files{n});
                fprintf('Deleting file %s \n',FileName)
                delete(FileName)
            end
        end
        UpdateFileTable(hObject,handles);
        guidata(hObject,handles)
    end
else 
    errordlg('No files selected!')
end


function DirectoryText_Callback(hObject, eventdata, handles)
% hObject    handle to DirectoryText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DirectoryText as text
%        str2double(get(hObject,'String')) returns contents of DirectoryText as a double
% Do not allow edit!
set(handles.DirectoryText,'String',handles.OutputDirectory)


% --- Executes during object creation, after setting all properties.
function DirectoryText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DirectoryText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

% --- Executes on button press in SelectDirectoriesButton.
function SelectDirectoriesButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectDirectoriesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
D=SelectDirectoryGui;       % Let user select those directories that (may) contain .imptime files to merge
% Now, for each directory, find .imptime files, and merge into the output
% directory.
for n = 1:length(D)
    ImptimeFiles=dir(fullfile(D{n},'*.imptime'));
    for m = 1:length(ImptimeFiles)
        MergeFile(handles.OutputDirectory,fullfile(D{n},ImptimeFiles(m).name));
    end
end
UpdateFileTable(hObject,handles);
guidata(hObject,handles)


%=========================================================================
% User Functions
function UpdateFileTable(hObject,handles)
D=dir(fullfile(handles.OutputDirectory,'*.imptime'));
FileNames={D.name};
for n=1:length(FileNames)
    load(fullfile(handles.OutputDirectory,FileNames{n}),'-mat');
    List{n,1}=FileNames{n};
    if length(Data) == 2
        List{n,2}=strcat('(',char(headpstr(Header(2))),')/',char(headpstr(Header(1))));
    else
        List{n,2}=strcat('*/',char(headpstr(Header(1))));
    end
end
if ~exist('List','var')
    List{1,1}='Empty Table';
end
set(handles.FileTable,'Data',List)
guidata(hObject,handles)


function MergeFile(OutputDirectory,FileName)
% Merge one file to the target directory. Check if a file with same name exists, and change if so.
% FileName includes full path

% First check if FileName includes Prefix AND number (if only Prefix, add
% the number '1')
InFileName=FileName;        % Original file name, needed later
[Path,Prefix,Number,Ext] = asplitfilename(FileName);
if isempty(Number)
    Number=1;
    FileName = strcat(Prefix,'1',Ext);
else
    FileName=strcat(Prefix,int2str(Number),Ext);
end
% Next, check if a file with same name exist in the output directory. If
% so, increment file number until a file number which does not already
% exist in the output directory is found. Then copy the file. If the
% filename does not exist, copy directly.
if exist(fullfile(OutputDirectory,FileName),'file') == 2
    % Increment file name; first find what current files exist in
    % terms of prefix and number. NOTE! It is assumed that all
    % files in the output directory conform to the Prefix<number>
    % convention, all with the same prefix. Best is to start with an empty
    % output directory and merge all your files there.
    Number=Number+1;
    FileName=strcat(Prefix,num2str(Number),Ext);
    while exist(fullfile(OutputDirectory,FileName),'file') == 2
        Number=Number+1;
        FileName=strcat(Prefix,num2str(Number),Ext);
    end
    OutFileName=fullfile(OutputDirectory,FileName);
    fprintf('Saving file with new name: %s\n',FileName)
    copyfile(InFileName,OutFileName);
else
    fprintf('Saving file %s\n',FileName)
    copyfile(InFileName,fullfile(OutputDirectory,FileName));
end


% --- Executes when selected cell(s) is changed in FileTable.
function FileTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to FileTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
handles.selection = eventdata.Indices(:,1);
guidata(hObject,handles);
