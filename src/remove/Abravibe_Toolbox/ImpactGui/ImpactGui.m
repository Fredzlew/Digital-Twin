function varargout = ImpactGui(varargin)
% IMPACTGUI ABRAVIBE impact testing GUI
%
%      h = ImpactGui opens the main GUI which controls the entire post
%      processing. The handle variable, h, is optional
%
%      This is a comprehensive, GUI-driven application that processes time
%      data from impact test recordings. The files with the time data
%      should have extension .imptime, and have the file format
%      documented in the ABRAVIBE User Manual. 
%
%      Several different impact testing process modes are supported, see
%      help for the command Imp2FrfGui. You can easily plug your own
%      processing mode into this file if you like.
%
%      The last settings used in the GUI are saved in a file in the abravibe 
%      directory called ImpactGuiSave.mat. If you want to revert to default 
%      settings for some reason, you can delete that file.
%
%      NOTE! Since ABRAVIBE is aimed at learning vibration analysis,
%      seemingly 'stupid' things, like selecting double impacts, are
%      allowed. 
%
%      Also have a look at the help documentation available directly from 
%      this GUI.
%
% See also: ABRA2IMP IMPACTGUIRESET

% Copyright (c) 2009-2015 by Anders Brandt
% 
% Version: 1.0 2014-07-15
% This file is part of ABRAVIBE Toolbox for NVA

% Last Modified by GUIDE v2.5 18-Dec-2015 11:51:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImpactGui_OpeningFcn, ...
                   'gui_OutputFcn',  @ImpactGui_OutputFcn, ...
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


% --- Executes just before ImpactGui is made visible.
function ImpactGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImpactGui (see VARARGIN)

% Choose default command line output for ImpactGui
handles.output = hObject;

%==============================================================================
% Some words about the impementation:
% Fields directly under the handles struct are used for global variables
% necessary only in this GUI (temporarily). A substruct field handles.ImpactGui 
% is used for all variables that are stored when you issue "save settings".
%==============================================================================
% Define fixed variables
% ONLY set fixed variables here...
handles.ImpactGui.BlocksizeAlfaList={'0.5K','1K','2K','4K','8K','16K','32K','64K','128K','256K'};
handles.ImpactGui.BlocksizeNumbers=1024*[0.5 1 2 4 8 16 32 64 128 256];
handles.ImpactGui.DefaultBlocksizeNumber=3;  % 2K default blocksise!
%==============================================================================

% Define default field values in handles struct
handles.SetupFile=which('ImpactGuiSave.mat');
if ~isempty(handles.SetupFile)
    load(handles.SetupFile)
    handles.ImpactGui=ImpactStruct;
else
    D=which('ImpactGui.m');
    [Path,Fn,Ext]=fileparts(D);    % Next clear file path from filename
    D=strcat(Fn,Ext);
    handles.ImpactGuiPath=Path;
    handles.SetupFile=fullfile(Path,'ImpactGuiSave');
    F=dir('*.imptime');
    if ~isempty(F)
        Fnames={F.name};
    else
        Fnames={''};
    end
    handles.ImpactGui.FunctionPath=Path;       % Points to ImpactGui directory
    handles.ImpactGui.InputDirectory=pwd;
    handles.ImpactGui.InputFileList=Fnames;
    handles.ImpactGui.OutputDirectory=fullfile(handles.ImpactGui.InputDirectory);
%     if exist(handles.ImpactGui.OutputDirectory,'dir') ~= 7
%         mkdir(handles.ImpactGui.OutputDirectory);
%     end
    handles.ImpactGui.OutputPrefix='abra';
    handles.ImpactGui.OutputStartNo=1;
    handles.ImpactGui.PretestFileName=1;
    handles.ImpactGui.PretestChannel=2;
    handles.ImpactGui.TrigLevel=5;
    handles.ImpactGui.Pretrigger=100;
    handles.ImpactGui.BlocksizeNumber=handles.ImpactGui.DefaultBlocksizeNumber;
    handles.ImpactGui.ForceWindowLength=100;
    handles.ImpactGui.ExpWindowPercent=100;
    handles.ImpactGui.ProcessMode='Manual';
    handles.ImpactGui.UseChannelNo=2;
end
    
% Fill GUI opening values
% Show only 39 last characters in input directory string
L=length(handles.ImpactGui.InputDirectory);
set(handles.InputDirectoryText,'string',['...' handles.ImpactGui.InputDirectory(L+1-min(L,39):L)]);
if strcmp(handles.ImpactGui.InputFileList{1},'')
    Fstring=sprintf('No files selected!',length(handles.ImpactGui.InputFileList));
else
    Fstring=sprintf('%i files selected!',length(handles.ImpactGui.InputFileList));
end
set(handles.InputNumberFilesText,'string',Fstring);
set(handles.OutputDirectoryEdit,'string',handles.ImpactGui.OutputDirectory);
set(handles.OutputPrefixEdit,'string',handles.ImpactGui.OutputPrefix);
set(handles.OutputStartNoEdit,'string',handles.ImpactGui.OutputStartNo);
set(handles.ContinueText,'String','When ready with pretest and settings click Continue to open up Gui for processing data!');
set(handles.PretestUseFilePopup,'Value',1);

% Fill popup menu for imptime files in input directory based on Input Files
% panel values
FillUseFilesPopoup(handles)
% Fill the Settings panel
set(handles.TrigLevelText,'string',num2str(handles.ImpactGui.TrigLevel));
set(handles.PretriggerText,'string',num2str(handles.ImpactGui.Pretrigger));
set(handles.BlocksizeText,'string',num2str(handles.ImpactGui.BlocksizeAlfaList{handles.ImpactGui.BlocksizeNumber}));
set(handles.ForceWindowText,'string',num2str(handles.ImpactGui.ForceWindowLength));
set(handles.ExpWindowText,'string',num2str(handles.ImpactGui.ExpWindowPercent));

% Define help menu and items
h1 = uimenu('Label', 'Help');
     uimenu(h1,'Label','Help','Callback','ImpHelp');
     uimenu(h1,'Label','About','Callback','ImpAbout');

guidata(hObject, handles);

% UIWAIT makes ImpactGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImpactGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SelectedImpactsListbox.
function SelectedImpactsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to SelectedImpactsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SelectedImpactsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SelectedImpactsListbox


% --- Executes during object creation, after setting all properties.
function SelectedImpactsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectedImpactsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ChangeOutputDirectoryPushbutton.
function ChangeOutputDirectoryPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeOutputDirectoryPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function OutputDirectoryEdit_Callback(hObject, eventdata, handles)
% hObject    handle to OutputDirectoryEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputDirectoryEdit as text
%        str2double(get(hObject,'String')) returns contents of OutputDirectoryEdit as a double
D=get(hObject,'String');
while exist(D,'dir') ~= 7
    D=uigetdir('.','Directory does not exist! Select and existing directory!');
end
set(hObject,'String',D);
handles.ImpactGui.OutputDirectory=D;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function OutputDirectoryEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputDirectoryEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SelectFilesButton.
function SelectFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Fn,Path]=uigetfile('*.imptime',...
    'Select files for analysis!',handles.ImpactGui.InputDirectory,'MultiSelect','on');
if iscell(Fn)
    handles.ImpactGui.InputDirectory=Path;
    L=length(handles.ImpactGui.InputDirectory);
    handles.ImpactGui.InputFileList=Fn;
    S=sprintf('%i files selected!',length(handles.ImpactGui.InputFileList));
    set(handles.InputDirectoryText,'string',['...' handles.ImpactGui.InputDirectory(L+1-min(L,39):L)]);
    set(handles.InputNumberFilesText,'String',S);
elseif isstr(Fn)
    handles.ImpactGui.InputDirectory=Path;
    L=length(handles.ImpactGui.InputDirectory);
    handles.ImpactGui.InputFileList={Fn};
    S=sprintf('%i files selected!',length(handles.ImpactGui.InputFileList));
    set(handles.InputDirectoryText,'string',['...' handles.ImpactGui.InputDirectory(L+1-min(L,39):L)]);
    set(handles.InputNumberFilesText,'String',S);
else
    ImpactContDlg('Title','No files selected! Please select files to proceed!');
end
FillUseFilesPopoup(handles)
guidata(hObject,handles);

% --- Executes on button press in CopyFilesButton.
function CopyFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to CopyFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CopyFilesGui('dummy',pwd,handles.ImpactGui.OutputDirectory);


function TrigLevelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TrigLevelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrigLevelEdit as text
%        str2double(get(hObject,'String')) returns contents of TrigLevelEdit as a double
n = str2double(get(hObject,'String'));
if isnumeric(n)
    handles.ImpactGui.TrigLevel=n;
    guidata(hObject,handles)
else
    ImpactContDlg('Title','Trig Level must be numeric! Not changing...')
    set(hObject,'String',num2str(handles.ImpactGui.TrigLevel));
end


% --- Executes during object creation, after setting all properties.
function TrigLevelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrigLevelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PretriggerEdit_Callback(hObject, eventdata, handles)
% hObject    handle to PretriggerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PretriggerEdit as text
%        str2double(get(hObject,'String')) returns contents of PretriggerEdit as a double
n = str2double(get(hObject,'String'));
if isnumeric(n)
    handles.ImpactGui.Pretrigger=n;
    guidata(hObject,handles)
else
    ImpactContDlg('Title','Trig Level must be numeric! Not changing...')
    set(hObject,'String',num2str(handles.ImpactGui.Pretrigger));
end


% --- Executes during object creation, after setting all properties.
function PretriggerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PretriggerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ForceWindowEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ForceWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ForceWindowEdit as text
%        str2double(get(hObject,'String')) returns contents of ForceWindowEdit as a double
n = str2double(get(hObject,'String'));
if isnumeric(n)
    if n > 0 & n <= 100
        handles.ImpactGui.ForceWindowLength=n;
        guidata(hObject,handles)
    else
        ImpactContDlg('Title','Force window length must be between 0 and 100%!')
        set(hObject,'String',num2str(handles.ImpactGui.ForceWindowLength));
    end
else
    ImpactContDlg('Title','Force window length must be numeric! Not changing...')
    set(hObject,'String',num2str(handles.ImpactGui.ForceWindowLength));
end



% --- Executes during object creation, after setting all properties.
function ForceWindowEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ForceWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ExpWindowEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ExpWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExpWindowEdit as text
%        str2double(get(hObject,'String')) returns contents of ExpWindowEdit as a double
n = str2double(get(hObject,'String'));
if isnumeric(n)
    if n > 0 & n <= 100
        handles.ImpactGui.ExpWindowPercent=n;
        guidata(hObject,handles)
    else
        ImpactContDlg('Title','Exponential window End must be between 0 and 100%!')
        set(hObject,'String',num2str(handles.ImpactGui.ExpWindowEndPercent));
    end
else
    ImpactContDlg('Title','Exponential window end must be numeric! Not changing...')
    set(hObject,'String',num2str(handles.ImpactGui.ExpWindowEndPercent));
end


% --- Executes during object creation, after setting all properties.
function ExpWindowEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExpWindowEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BlocksizePopup.
function BlocksizePopup_Callback(hObject, eventdata, handles)
% hObject    handle to BlocksizePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BlocksizePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BlocksizePopup
contents=cellstr(get(hObject,'String'));
handles.ImpactGui.BlocksizeNumber=get(hObject,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function BlocksizePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BlocksizePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PretestStartButton.
function PretestStartButton_Callback(hObject, eventdata, handles)
% hObject    handle to PretestStartButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.ImpactGui = ImpactSettings('String',handles.ImpactGui);
% First check if the popup menu for file has been used to select a file
n=get(handles.PretestUseFilePopup,'Value');
if n ~= 1 
    m=get(handles.PretestUseChannelPopup,'Value');
    if m > 1
        A=ImpactSettings(handles.ImpactGui);
        % Set GUI fields based on new settings if not empty fields from
        % ImpactSettings
        if ~isempty(A)
            handles.ImpactGui=A;
            set(handles.TrigLevelText,'string',num2str(handles.ImpactGui.TrigLevel));
            set(handles.PretriggerText,'string',num2str(handles.ImpactGui.Pretrigger));
            set(handles.BlocksizeText,'string',num2str(handles.ImpactGui.BlocksizeAlfaList{handles.ImpactGui.BlocksizeNumber}));
            set(handles.ForceWindowText,'string',num2str(handles.ImpactGui.ForceWindowLength));
            set(handles.ExpWindowText,'string',num2str(handles.ImpactGui.ExpWindowPercent));
            UpdateSetupFile(handles);
            guidata(hObject,handles);
        end
    else
        ImpactContDlg('Title','Channel to use is not selected!')
    end
else
        ImpactContDlg('Title','Select file and channel to use!') 
end

% --- Executes on button press in AutomaticProcessAllButton.
function AutomaticProcessAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to AutomaticProcessAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LoadSettingsButton.
function LoadSettingsButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadSettingsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ImpactGui.OutputDirectory
[FileName,PathName]=uigetfile(fullfile(handles.ImpactGui.OutputDirectory,...
    '*.ImpGuiSetup'),'Select setup file to open!');
if FileName ~= 0
    load(fullfile(PathName,FileName),'-mat');
    handles.ImpactGui=ImpactStruct;
end
% Update GUI with new settings
% set(handles.InputDirectoryEdit,'string',handles.ImpactGui.InputDirectory);
% set(handles.InputPrefixEdit,'string',handles.ImpactGui.InputPrefix);
% set(handles.InputStartNoEdit,'string',handles.ImpactGui.InputStartNo);
% set(handles.InputLastNoEdit,'string',handles.ImpactGui.InputLastNo);
set(handles.OutputDirectoryEdit,'string',handles.ImpactGui.OutputDirectory);
set(handles.OutputPrefixEdit,'string',handles.ImpactGui.OutputPrefix);
set(handles.OutputStartNoEdit,'string',handles.ImpactGui.OutputStartNo);
set(handles.ContinueText,'String','When ready with pretest and settings click Continue to open up Gui for processing data!');
set(handles.PretestUseFilePopup,'Value',1);
% Fill popup menu for imptime files in input directory based on Input Files
% panel values
FillUseFilesPopoup(handles)
set(handles.TrigLevelText,'string',num2str(handles.ImpactGui.TrigLevel));
set(handles.PretriggerText,'string',num2str(handles.ImpactGui.Pretrigger));
set(handles.BlocksizeText,'string',num2str(handles.ImpactGui.BlocksizeAlfaList{handles.ImpactGui.BlocksizeNumber}));
set(handles.ForceWindowText,'string',num2str(handles.ImpactGui.ForceWindowLength));
set(handles.ExpWindowText,'string',num2str(handles.ImpactGui.ExpWindowPercent));
% Save results in backup file
UpdateSetupFile(handles);

guidata(hObject,handles);

% --- Executes on button press in SaveSettingsButton.
function SaveSettingsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSettingsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName]=uiputfile(fullfile(handles.ImpactGui.OutputDirectory,...
    '*.ImpGuiSetup'),'Select file to save setup to!');
if FileName ~= 0
    FileName=fullfile(PathName,FileName);
    ImpactStruct=handles.ImpactGui;
    save(FileName,'ImpactStruct');
end


function OutputStartNoEdit_Callback(hObject, eventdata, handles)
% hObject    handle to OutputStartNoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputStartNoEdit as text
%        str2double(get(hObject,'String')) returns contents of OutputStartNoEdit as a double
n=str2double(get(hObject,'String'));
if ~isnan(n)
    handles.ImpactGui.OutputStartNo=n;
    guidata(hObject,handles)
else
    ImpactContDlg('Title','StartNo must not be empty! Not changing...') 
    set(hObject,'String',num2str(handles.ImpactGui.InputStartNo));
end


% --- Executes during object creation, after setting all properties.
function OutputStartNoEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputStartNoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OutputPrefixEdit_Callback(hObject, eventdata, handles)
% hObject    handle to OutputPrefixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OutputPrefixEdit as text
%        str2double(get(hObject,'String')) returns contents of OutputPrefixEdit as a double
s=get(hObject,'String');
if ~isempty(s)
    handles.ImpactGui.OutputPrefix=s;
    guidata(hObject,handles)
else
    ImpactContDlg('Title','Prefix must not be empty! Not changing...')
    set(hObject,'String',handles.ImpactGui.OutputPrefix);
end


% --- Executes during object creation, after setting all properties.
function OutputPrefixEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputPrefixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ChangeOutputDirectoryButton.
function ChangeOutputDirectoryButton_Callback(hObject, eventdata, handles)
% hObject    handle to ChangeOutputDirectoryButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
F=uigetdir(handles.ImpactGui.OutputDirectory,'Set input directory!');
if F ~= 0
    handles.ImpactGui.OutputDirectory=F;
    set(handles.OutputDirectoryEdit,'String',handles.ImpactGui.OutputDirectory);
    UpdateSetupFile(handles);
    guidata(hObject,handles)
end


% --- Executes on selection change in PretestUseFilePopup.
function PretestUseFilePopup_Callback(hObject, eventdata, handles)
% hObject    handle to PretestUseFilePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PretestUseFilePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PretestUseFilePopup
contents=cellstr(get(hObject,'String'));
FileNo=get(hObject,'Value');
if FileNo ~= 1
    FileName=fullfile(handles.ImpactGui.InputDirectory,handles.ImpactGui.InputFileList{FileNo-1});
    load(FileName,'-mat')
    NumberChannels=length(Data);
    ChannelList={'Use Channel'};
    for n=2:NumberChannels 
        ChannelList=[ChannelList;int2str(n)];
    end
%     set(handles.PretestUseChannelPopup,'String',[2:NumberChannels]);
    set(handles.PretestUseChannelPopup,'String',ChannelList);
    set(handles.PretestUseChannelPopup,'Value',1);
    handles.ImpactGui.PretestFileName=FileName;
    set(handles.PretestUseChannelPopup,'Visible','on')
    guidata(hObject,handles)
else
    set(handles.PretestUseChannelPopup,'Visible','off')
end

% --- Executes during object creation, after setting all properties.
function PretestUseFilePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PretestUseFilePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PretestUseChannelPopup.
function PretestUseChannelPopup_Callback(hObject, eventdata, handles)
% hObject    handle to PretestUseChannelPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PretestUseChannelPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PretestUseChannelPopup
contents=cellstr(get(hObject,'String'));
FileNo=get(hObject,'Value');
if FileNo ~= 1
    handles.ImpactGui.UseChannelNo=FileNo;
    guidata(hObject,handles)
end

% --- Executes during object creation, after setting all properties.
function PretestUseChannelPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PretestUseChannelPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InputPrefixEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InputPrefixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputPrefixEdit as text
%        str2double(get(hObject,'String')) returns contents of InputPrefixEdit as a double
s=get(hObject,'String');
if ~isempty(s)
    handles.ImpactGui.InputPrefix=s;
    FillUseFilesPopoup(handles)
    guidata(hObject,handles)
else
    ImpactContDlg('Title','Prefix must not be empty! Not changing...')
    set(hObject,'String',handles.ImpactGui.InputPrefix);
end

% --- Executes during object creation, after setting all properties.
function InputPrefixEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputPrefixEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InputStartNoEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InputStartNoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputStartNoEdit as text
%        str2double(get(hObject,'String')) returns contents of InputStartNoEdit as a double
n=str2double(get(hObject,'String'));
if ~isnan(n)
    if n > handles.ImpactGui.InputLastNo
        ImpactContDlg('Title','StartNo must be <= LastNo! Setting = LastNo')
        set(handles.InputStartNoEdit,'String',handles.ImpactGui.InputLastNo);
        handles.ImpactGui.InputStartNo=handles.ImpactGui.InputLastNo;
    else
        handles.ImpactGui.InputStartNo=n;
    end
    FillUseFilesPopoup(handles)
    guidata(hObject,handles)
else
    ImpactContDlg('Title','StartNo must not be empty! Not changing...') 
    set(hObject,'String',num2str(handles.ImpactGui.InputStartNo));
end


% --- Executes during object creation, after setting all properties.
function InputStartNoEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputStartNoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function InputLastNoEdit_Callback(hObject, eventdata, handles)
% hObject    handle to InputLastNoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputLastNoEdit as text
%        str2double(get(hObject,'String')) returns contents of InputLastNoEdit as a double
n=str2double(get(hObject,'String'));
if ~isnan(n)
    if n < handles.ImpactGui.InputStartNo
        ImpactContDlg('Title','LastNo must be <= StartNo! Setting = StartNo')
        set(handles.InputLastNoEdit,'String',handles.ImpactGui.InputStartNo);
        handles.ImpactGui.InputLastNo=handles.ImpactGui.InputStartNo;
    else
        handles.ImpactGui.InputLastNo=n;
    end
    FillUseFilesPopoup(handles)
    guidata(hObject,handles)
else
    ImpactContDlg('Title','LastNo must not be empty! Not changing...')
    set(hObject,'String',num2str(handles.ImpactGui.InputLastNo));
end


% --- Executes during object creation, after setting all properties.
function InputLastNoEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputLastNoEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ManualProcessButton.
function ManualProcessButton_Callback(hObject, eventdata, handles)
% hObject    handle to ManualProcessButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ManualNextButton.
function ManualNextButton_Callback(hObject, eventdata, handles)
% hObject    handle to ManualNextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in ContinueButton.
function ContinueButton_Callback(hObject, eventdata, handles)
% hObject    handle to ContinueButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ImpactProcessingGui(handles);

% --- Executes on selection change in PlotChannelsListbox.
function PlotChannelsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to PlotChannelsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlotChannelsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotChannelsListbox


% --- Executes during object creation, after setting all properties.
function PlotChannelsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotChannelsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in UseImpactsListbox.
function UseImpactsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to UseImpactsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UseImpactsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UseImpactsListbox


% --- Executes during object creation, after setting all properties.
function UseImpactsListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UseImpactsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%=====================================================================
% User defined functions
%=====================================================================

function FillUseFilesPopoup(handles)
% After an update in the input files settings, the popup menu with all
% imptime files in the current input directory with Prefix and StartNo,
% LastNo needs to be filled. If no files found, that is written to the
% popup menu
% Fill popup menu for imptime files in input directory

FileList={'Select File'};
FileList=[FileList handles.ImpactGui.InputFileList];
if length(FileList) == 1
    FileList={'No Files!'};
    set(handles.PretestUseFilePopup,'Value',1);
end
set(handles.PretestUseFilePopup,'String',FileList)
set(handles.PretestUseChannelPopup,'Visible','off')
ImpactStruct=handles.ImpactGui;
UpdateSetupFile(handles);


function UpdateSetupFile(handles)
% Save current data to setup file
ImpactStruct=handles.ImpactGui;
save(handles.SetupFile,'ImpactStruct');


% --- Executes on button press in BrowseFilesButton.
function BrowseFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
abrabrowse
