function varargout = ImpactSettings(varargin)
% IMPACTSETTINGS MATLAB code for ImpactSettings.fig
%      IMPACTSETTINGS, by itself, creates a new IMPACTSETTINGS or raises the existing
%      singleton*.
%
%      H = IMPACTSETTINGS returns the handle to a new IMPACTSETTINGS or the handle to
%      the existing singleton*.
%
%      IMPACTSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPACTSETTINGS.M with the given input arguments.
%
%      IMPACTSETTINGS('Property','Value',...) creates a new IMPACTSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImpactSettings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImpactSettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2009-2014 by Anders Brandt
% Email: abra@iti.sdu.dk
% Version: 1.0 2014-07-15
% This file is part of ABRAVIBE Toolbox for NVA


% Last Modified by GUIDE v2.5 26-Nov-2015 17:57:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ImpactSettings_OpeningFcn, ...
    'gui_OutputFcn',  @ImpactSettings_OutputFcn, ...
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


% --- Executes just before ImpactSettings is made visible.
function ImpactSettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImpactSettings (see VARARGIN)

% Define help menu and items
h1 = uimenu('Label', 'Help','Callback','ImpSettingsHelp');
%      uimenu(h1,'Label','Help','Callback','ImpSettingsHelp');

% Load input variables from file ImpactSettings
handles.ImpactGui=varargin{1};

% Set GUI fields by input values from the ImpactGui struct
set(handles.TrigLevelEdit,'String',handles.ImpactGui.TrigLevel);
set(handles.PretriggerEdit,'String',handles.ImpactGui.Pretrigger);
set(handles.BlocksizePopup,'Value',handles.ImpactGui.BlocksizeNumber);
set(handles.ForceWindowEdit,'String',handles.ImpactGui.ForceWindowLength);
set(handles.ExpWindowEdit,'String',handles.ImpactGui.ExpWindowPercent);

% Choose default command line output for ImpactSettings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImpactSettings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImpactSettings_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% Start by closing down all plot windows that are open
if isfield(handles,'h1') & ishghandle(handles.h1)
    delete(handles.h1)
    rmfield(handles,'h1')
end
if isfield(handles,'h2') & ishghandle(handles.h2)
    delete(handles.h2)
    rmfield(handles,'h2')
end
if isfield(handles,'h3') & ishghandle(handles.h3)
    delete(handles.h3)
    rmfield(handles,'h3')
end
if isfield(handles,'h4') & ishghandle(handles.h4)
    delete(handles.h4)
    rmfield(handles,'h4')
end
if isfield(handles,'ImpactGui')
    handles.output=handles.ImpactGui;
    varargout{1} = handles.output;
    delete(handles.figure1);
else
    varargout{1} ='';
end

% --- Executes on button press in Export.
function Export_Callback(hObject, eventdata, handles)
% hObject    handle to Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;


% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ImpactGui='';
guidata(hObject,handles)
uiresume;


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


% --- Executes on button press in WindowsPlotFRFButton.
function WindowsPlotFRFButton_Callback(hObject, eventdata, handles)
% hObject    handle to WindowsPlotFRFButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'TrigIdx')
    N=handles.ImpactGui.BlocksizeNumbers(handles.ImpactGui.BlocksizeNumber);
    [H,f, C, Tff] = imp2frf2(handles.x,handles.y,handles.fs,N,handles.TrigIdx,...
        handles.DIdx,handles.ImpactGui.ForceWindowLength,...
        handles.ImpactGui.ExpWindowPercent,0);
    if ~isfield(handles,'h3')
        handles.h3=figure;
    elseif ishghandle(handles.h3)
        figure(handles.h3);
    else
        handles.h3=figure;
    end
    subplot(3,1,1)
    semilogy(f,abs(H))
    ylabel('Mag. FRF')
    title('Zoom in and investigate the lowest resonance peak(s)')
    axis tight
    A=axis;
    subplot(3,1,2)
    plot(f,C)
    ylabel('Coherence')
    xlim(A(1:2))
    ylim([0 1.02])
    subplot(3,1,3)
    semilogy(f,Tff)
    xlim(A(1:2))
    ylabel('Force Spectrum')
    xlabel('Frequency [Hz]')
    guidata(hObject,handles)
else
    ImpactContDlg('Title','You have to select impacts first!')
end


% --- Executes on selection change in BlocksizePopup.
function BlocksizePopup_Callback(hObject, eventdata, handles)
% hObject    handle to BlocksizePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BlocksizePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BlocksizePopup
contents=cellstr(get(hObject,'String'));
OldBlocksizeNumber=handles.ImpactGui.BlocksizeNumber;
handles.ImpactGui.BlocksizeNumber=get(hObject,'Value');
% If selected immediately upon opening GUI, there are no triggers, so only
% update the blocksize. Otherwise select new triggers, and compare with
% previous triggers
if isfield(handles,'x')     % Indicates triggers are already selected
    % Check if new blocksize affects the triggers or double impact detection
    % First get new triggers
    TrigLevel=handles.ImpactGui.TrigLevel;
    Pretrigger=handles.ImpactGui.Pretrigger;
    N=handles.ImpactGui.BlocksizeNumbers(handles.ImpactGui.BlocksizeNumber);
    [TrigIdx,DIdx,AIdx]=GuiImptrig(handles.x,N,TrigLevel,Pretrigger);
    % Check if new trigger occations are different from before
    Idx=ismember(handles.TrigIdx,TrigIdx);
    if ~isempty(find(Idx == 0))   % if trigger occasions have changed
        if isfield(handles,'TrigIdx')
            handles=rmfield(handles,'TrigIdx');
            set(handles.NumberImpactsText,'String','Triggering affected! Select new impacts!')
        end
    end
    % Check that no blocks are overlapping with the new blocksize, unless
    % trigger index already deleted
    if isfield(handles,'TrigIdx')
        [status,TrigIdx]=ImpactChecks(handles.x,handles.TrigIdx,N);
        if ~isempty(TrigIdx)
            handles.TrigIdx = TrigIdx;
            set(handles.NumberImpactsText,'String',[int2str(length(TrigIdx)) ' impacts selected'])
        else
            set(handles.NumberImpactsText,'String','Select New Impacts!')
            %         handles.ImpactGui.BlocksizeNumber=OldBlocksizeNumber;
        end
    end
end    
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


% --- Executes on button press in BlocksizeButton.
function BlocksizeButton_Callback(hObject, eventdata, handles)
% hObject    handle to BlocksizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Plot FRFs with three blocksizes
if isfield(handles,'TrigIdx')
    N=handles.ImpactGui.BlocksizeNumbers(handles.ImpactGui.BlocksizeNumber);
    N=N*[0.5 1 2];
    % Make sure the larger blocksize does not exhaust data
    [status,TrigIdx]=ImpactChecks(handles.x,handles.TrigIdx,N(3));
    if isempty(TrigIdx)
        NumberBlocksizes=2;
    else
        NumberBlocksizes=3;
    end
    for n = 1:NumberBlocksizes
        % Compute FRFs. NOTE! This only works if you use no windows...
        [H{n},f{n}] = imp2frf2(handles.x,handles.y,handles.fs,N(n),handles.TrigIdx,...
            handles.DIdx,100,100,0);
        if n==1
            L={int2str(N(n))};
        else
            L=[L; int2str(N(n))];
        end
    end
    if ~isfield(handles,'h2')
        handles.h2=figure;
    elseif ishghandle(handles.h2)
        figure(handles.h2);
    else
        handles.h2=figure;
    end
    semilogy(f{1},abs(H{1}),f{2},abs(H{2}),f{3},abs(H{3}))
    legend(L)
    xlabel('Frequency [Hz]')
    ylabel('Mag. FRF')
    title('Zoom in and investigate the lowest resonance peak(s)')
    guidata(hObject,handles)
else
    ImpactContDlg('Title','You have to select impacts first!')
end



function TrigLevelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TrigLevelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrigLevelEdit as text
%        str2double(get(hObject,'String')) returns contents of TrigLevelEdit as a double
n = str2double(get(hObject,'String'));
if isnumeric(n)
    handles.ImpactGui.TrigLevel=n;
    % Reset trigger index if it exists (i.e. if impacts have been selected)
    if isfield(handles,'TrigIdx')
        handles=rmfield(handles,'TrigIdx');
        set(handles.NumberImpactsText,'String','No Impacts Selected!')
    end
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
    % Reset trigger index if it exists (i.e. if impacts have been selected)
    if isfield(handles,'TrigIdx')
        handles=rmfield(handles,'TrigIdx');
        set(handles.NumberImpactsText,'String','No Impacts Selected!')
    end
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


% --- Executes on button press in TriggeringButton.
function TriggeringButton_Callback(hObject, eventdata, handles)
% hObject    handle to TriggeringButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'TrigIdx')
    % Plot force and response signals in blue. Mark the "good" length N blocks
    % with no double impacts in green, and double impact blocks in red
    N=handles.ImpactGui.BlocksizeNumbers(handles.ImpactGui.BlocksizeNumber);
    t=makexaxis(handles.x,1/handles.fs);
    if ~isfield(handles,'h1') 
        handles.h1=figure;
    elseif ishghandle(handles.h1) 
        figure(handles.h1);
    else
        handles.h1=figure;
    end
    subplot(2,1,1)
    TrigLevelScaled=max(handles.x)*handles.ImpactGui.TrigLevel/100;
    plot(t,handles.x,[t(1) t(end)],[TrigLevelScaled TrigLevelScaled])
    title('Zoom in and investigate impacts and response! Then close window.')
    hold on
    for n=1:length(handles.TrigIdx)
        idx=handles.TrigIdx(n):handles.TrigIdx(n)+N-1;
        plot(t(idx),handles.x(idx),'g')
    end
    if handles.DIdx ~= 0
        for n = 1:length(handles.DIdx)
            idx=handles.AIdx(handles.DIdx(n)):handles.AIdx(handles.DIdx(n))+N-1;
            plot(t(idx),handles.x(idx),'r')
        end
    end
    hold off
    subplot(2,1,2)
    plot(t,handles.y)
    hold on
    for n=1:length(handles.TrigIdx)
        idx=handles.TrigIdx(n):handles.TrigIdx(n)+N-1;
        plot(t(idx),handles.y(idx),'g')
    end
    if handles.DIdx ~= 0
        for n = 1:length(handles.DIdx)
            idx=handles.AIdx(handles.DIdx(n)):handles.AIdx(handles.DIdx(n))+N-1;
            plot(t(idx),handles.y(idx),'r')
        end
    end
    guidata(hObject,handles)
else
    ImpactContDlg('Title','You have to select impacts first!')
end

% --- Executes on button press in SelectImpactsButton.
function SelectImpactsButton_Callback(hObject, eventdata, handles)
% hObject    handle to SelectImpactsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load(handles.ImpactGui.PretestFileName,'-mat')
handles.x=Data{1};
% Check if impacts are negative, then "flip" positive, and update Dir field
handles.x=cell2mat(Data(1));
[M,I]=max(abs(handles.x));
if handles.x(I) < 0
    handles.x=-handles.x;
    Dir=dir2nbr(Header(1).Dir);
    Header.Dir=nbr2dir(-Dir);
end
handles.y=Data{handles.ImpactGui.UseChannelNo};
handles.fs=1/Header(1).xIncrement;
N=handles.ImpactGui.BlocksizeNumbers(handles.ImpactGui.BlocksizeNumber);
TrigLevel=handles.ImpactGui.TrigLevel;
Pretrigger=handles.ImpactGui.Pretrigger;
[TrigIdx,DIdx,AIdx]=GuiImptrig2(handles.x,handles.y,handles.fs,N,TrigLevel,Pretrigger);
handles.TrigIdx=TrigIdx;
handles.DIdx=DIdx;
handles.AIdx=AIdx;
set(handles.NumberImpactsText,'String',[int2str(length(TrigIdx)) ' impacts selected'])
% Check if there are double impacts among the selected impacts. In that
% case issue a warning. Thus, this is 'allowed', because you may want to see
% what the result is.
if DIdx ~= 0
    if ~isempty(intersect(TrigIdx,AIdx(DIdx)))
        ImpactContDlg('Title','WARNING! You have selected blocks with double impacts!')
    end
end
guidata(hObject,handles);


% --- Executes on button press in WindowsPlotTimeButton.
function WindowsPlotTimeButton_Callback(hObject, eventdata, handles)
% hObject    handle to WindowsPlotTimeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'TrigIdx')
    if ~isfield(handles,'h4')
        handles.h4=figure;
    elseif ishghandle(handles.h4)
        figure(handles.h4);
    else
        handles.h4=figure;
    end
    N=handles.ImpactGui.BlocksizeNumbers(handles.ImpactGui.BlocksizeNumber);
    idx=handles.TrigIdx(1):handles.TrigIdx(1)+N-1;
    t=makexaxis(handles.x(idx),1/handles.fs);
    wf=aforcew(N,handles.ImpactGui.ForceWindowLength);
    we=aexpw(N,handles.ImpactGui.ExpWindowPercent);
    subplot(3,1,1)
    plot(t,handles.x(idx),t,wf*max(handles.x(idx))*1.02);
    ylabel('Force, First Impact')
    subplot(3,1,2)
    plot(t,handles.y(idx))
    ylabel('Response Signal')
    Ma=max(abs(handles.y(idx)));
    ylim([-Ma Ma])
    subplot(3,1,3)
    plot(t,handles.y(idx).*we)
    ylabel('Windowed Response')
    xlabel('Time [s]')
    ylim([-Ma Ma])
    guidata(hObject,handles)
else
    ImpactContDlg('Title','You have to select impacts first!')
end


% --- Executes during object creation, after setting all properties.
function NumberImpactsText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumberImpactsText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
