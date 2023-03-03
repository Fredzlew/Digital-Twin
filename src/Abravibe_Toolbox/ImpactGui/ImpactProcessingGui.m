function varargout = ImpactProcessingGui(varargin)
% IMPACTPROCESSINGGUI MATLAB code for ImpactProcessingGui.fig
%      IMPACTPROCESSINGGUI, by itself, creates a new IMPACTPROCESSINGGUI or raises the existing
%      singleton*.
%
%      H = IMPACTPROCESSINGGUI returns the handle to a new IMPACTPROCESSINGGUI or the handle to
%      the existing singleton*.
%
%      IMPACTPROCESSINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPACTPROCESSINGGUI.M with the given input arguments.
%
%      IMPACTPROCESSINGGUI('Property','Value',...) creates a new IMPACTPROCESSINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImpactProcessingGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImpactProcessingGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2009-2014 by Anders Brandt
% Email: abra@iti.sdu.dk
% Version: 1.0 2014-07-15
%          1.1 2018-05-12 Fixed bug in active impacts after first save
% This file is part of ABRAVIBE Toolbox for NVA


% Last Modified by GUIDE v2.5 12-Jul-2014 15:05:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ImpactProcessingGui_OpeningFcn, ...
    'gui_OutputFcn',  @ImpactProcessingGui_OutputFcn, ...
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


% --- Executes just before ImpactProcessingGui is made visible.
function ImpactProcessingGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImpactProcessingGui (see VARARGIN)

% Keep handles to Main gui to be able to update file numbers in that GUI as
% files are processed
MainHandles=varargin{1};
% set(MainHandles.InputStartNoEdit,'String','3');
handles.ImpactGui=MainHandles.ImpactGui;
handles.ProcessType='FREQUENCYMANUAL';

% Load first data file based on the input settings and fill the Acc.
% Channels listbox
FileName=fullfile(handles.ImpactGui.InputDirectory,handles.ImpactGui.InputFileList{1});
if exist(FileName,'file') == 2
    load(FileName,'-mat')
    handles.Data=Data;
    handles.Header=Header;
else
    ImpactContDlg('Title','Input file does not exist!')
end
handles.CurrentFile=FileName;
NumberChannels=length(Data);    % Number of channels including force channel
PlotForArray=cell(NumberChannels,1);
PlotForArray{1}=['Channel(s)'];
PlotForArray{2}=[num2str(2,'%4i')];
for n=3:NumberChannels
    PlotForArray{n}=[num2str(n,'%4i')];
end
set(handles.PlotChannelsListbox,'String',PlotForArray)
handles.UseChannels=2;
set(handles.PlotChannelsListbox,'Value',handles.UseChannels)

% Now fill Use Impacts table
% For this we first need data
handles.x=cell2mat(Data(1));
handles.y=cell2mat(Data(2:end));
handles.fs=1/Header(1).xIncrement;
handles.N=handles.ImpactGui.BlocksizeNumbers(handles.ImpactGui.BlocksizeNumber);
% Now find trigger events, including double impacts
[handles.TrigIdx,handles.DIdx,handles.AIdx]=GuiImptrig(handles.x,handles.N,...
    handles.ImpactGui.TrigLevel,handles.ImpactGui.Pretrigger);
% Fill the listbox with as many impacts as were found
C=cell(length(handles.TrigIdx),1);
for n = 1:length(handles.TrigIdx)
    C{n} = num2str(n,'%4i');
end
set(handles.UseImpactsListbox,'String',C)
% Make only those not being double impacts as active. This allows the user
% to select double impact blocks, but they are not used by default
handles.UseImpacts=setdiff(1:length(handles.TrigIdx),handles.DIdx);
set(handles.UseImpactsListbox,'Value',handles.UseImpacts)

% Fill the Current file popup
FillUseFilesPopoup(handles);


% Next, plot FRF and coherence for default acc. channel no. 2 and using
% default (all) impacts
ScreenSize=get(0,'ScreenSize');
% Plot time plot
handles.h1=figure('Position',[ScreenSize(3)/2 ScreenSize(4)*2/3-15 ScreenSize(3)/2*.8 ScreenSize(4)/3*.8]);
UpdateTimePlots(handles);
% Plot mag. FRF and Coherence
handles.h2=figure('Position',[ScreenSize(3)/2 ScreenSize(2)+40 ScreenSize(3)/2*.8 ScreenSize(4)*2/3*.8]);
handles=UpdateFreqPlots(handles);

% Choose default command line output for ImpactProcessingGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImpactProcessingGui wait for user response (see UIRESUME)uiwait(handles.figure1);
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImpactProcessingGui_OutputFcn(~, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
if isfield(handles,'h1') & ishghandle(handles.h1)
    close(handles.h1)
end
if isfield(handles,'h2') & ishghandle(handles.h2)
    close(handles.h2)
end

% --- Executes on button press in SaveButton.
function SaveButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Save FRF, Coherence, and force spectrum for each response channel
% This means the force spectrum is repeated; this makes it easier to figure
% out which force spectrum belongs to which FRF
if strcmp(upper(handles.ProcessType),'FREQUENCYMANUAL') | ...
        strcmp(upper(handles.ProcessType),'TIMESYNCHRONOUSMANUAL')
    handles = SaveFiles(handles);
    % Update current file popup to next file, if not end of list!
    % If list is up, close Post Processing window
    FileNames=get(handles.CurrentFilePopup,'String');
    FileNo=get(handles.CurrentFilePopup,'Value');
    if FileNo < length(FileNames)
        FileNo=FileNo+1;
        set(handles.CurrentFilePopup,'Value',FileNo);
        FileName=fullfile(handles.ImpactGui.InputDirectory,FileNames{FileNo});
        load(FileName,'-mat');
        handles.CurrentFile=FileName;
        % Check if impacts are negative, then "flip" positive
        handles.x=cell2mat(Data(1));
        [M,I]=max(abs(handles.x));
        if handles.x(I) < 0
            handles.x=-handles.x;
            Dir=dir2nbr(Header(1).Dir);
            Header.Dir=nbr2dir(-Dir);
        end
        handles.y=cell2mat(Data(2:end));
        handles.Header=Header;
        [handles.TrigIdx,handles.DIdx,handles.AIdx]=GuiImptrig(handles.x,handles.N,...
            handles.ImpactGui.TrigLevel,handles.ImpactGui.Pretrigger);
        % Fill the listbox with as many impacts as were found
        C=cell(length(handles.TrigIdx),1);
        for n = 1:length(handles.TrigIdx)
            C{n} = num2str(n,'%4i');
        end
        set(handles.UseImpactsListbox,'String',C)
        %         handles.UseImpacts=[1:length(handles.TrigIdx)];
%         set(handles.UseImpactsListbox,'String',{handles.UseImpacts'});
        handles.UseImpacts=setdiff(1:length(handles.TrigIdx),handles.DIdx);
        set(handles.UseImpactsListbox,'Value',handles.UseImpacts);

        UpdateTimePlots(handles);
        handles=UpdateFreqPlots(handles);
        guidata(hObject,handles)
    else
        if isfield(handles,'h1') & ishghandle(handles.h1)
            close(handles.h1)
        end
        if isfield(handles,'h2') & ishghandle(handles.h2)
            close(handles.h2)
        end
        uiresume;
        delete(handles.figure1);
    end
%     guidata(hObject,handles)
else % if automatic processing mode: process and save results for all files
     % in the input list
    % Reread the current file in the UseFiles popup, and loop through all
    % files, process, and save
    FileNames=get(handles.CurrentFilePopup,'String');
    FileNo=get(handles.CurrentFilePopup,'Value');
    h=waitbar(0,'Processing .imptime files...');
    for n = FileNo:length(FileNames)
        waitbar(n/length(FileNames),h)
        figure(h)
        FileName=fullfile(handles.ImpactGui.InputDirectory,FileNames{n});
        load(FileName,'-mat')
        handles.CurrentFile=FileName;
        % Check if impacts are negative, then "flip" positive
        handles.x=cell2mat(Data(1));
        [M,I]=max(abs(handles.x));
        if handles.x(I) < 0
            handles.x=-handles.x;
            Dir=dir2nbr(Header(1).Dir);
            Header.Dir=nbr2dir(-Dir);
        end
        handles.y=cell2mat(Data(2:end));
        handles.Header=Header;
        % Compute new triggers
        [handles.TrigIdx,handles.DIdx,handles.AIdx]=GuiImptrig(handles.x,handles.N,...
            handles.ImpactGui.TrigLevel,handles.ImpactGui.Pretrigger);
        % Compute FRFs, Coherences, and force spectrum
        [handles.H,handles.f, handles.C, handles.Tff,handles.UseImpacts] = Imp2FrfGui(handles.x,...
            handles.y,handles.fs,handles.N,handles.TrigIdx,...
            handles.DIdx,handles.ImpactGui.ForceWindowLength,...
            handles.ImpactGui.ExpWindowPercent,handles.ProcessType,...
            handles.flo,handles.fhi);
        % Save files
        handles = SaveFiles(handles);
    end % for 
    % Close all windows
    if isfield(handles,'h1') & ishghandle(handles.h1)
        close(handles.h1)
    end
    if isfield(handles,'h2') & ishghandle(handles.h2)
        close(handles.h2)
    end
    uiresume;
    delete(handles.figure1);
    close(h)
end

% --- Executes on selection change in ProcessModePopup.
function ProcessModePopup_Callback(hObject, eventdata, handles)
% hObject    handle to ProcessModePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ProcessModePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ProcessModePopup
n=get(hObject,'Value');
if n == 1  % Freq. manual mode
    if isfield(handles,'flo')
        handles=rmfield(handles,'flo');
    end
    handles.ProcessType='FrequencyManual';
    set(handles.SaveButton,'String','Save and Next')
    handles=UpdateFreqPlots(handles);
    set(handles.PlotChannelsListbox,'Visible','on')
    set(handles.UseImpactsListbox,'Visible','on')
    set(handles.CurrentFilePopup,'Visible','on')
elseif n == 2 % Time Synchr. manual mode
    if isfield(handles,'flo')
        handles=rmfield(handles,'flo');
    end
    handles.ProcessType='TimeSynchronousManual';
    set(handles.SaveButton,'String','Save and Next')
    handles=UpdateFreqPlots(handles);
    set(handles.PlotChannelsListbox,'Visible','on')
    set(handles.UseImpactsListbox,'Visible','on')
    set(handles.CurrentFilePopup,'Visible','on')
elseif n == 3 % Best two optimization mode
    figure(handles.h2)
    subplot(4,1,[1 2])
    title('Select lower and upper frequency range for optimization!')
    h=helpdlg('Choose upper and lower frequencies in the FRF/Coherence plot! Then select OK!');
    figure(handles.h2)
    [xx,dum]=ginput(2);
    waitfor(h)
    handles.flo=xx(1);
    handles.fhi=xx(2);
    handles.ProcessType='BestTwo';
    set(handles.SaveButton,'String','Process All')
    set(handles.PlotChannelsListbox,'Visible','off')
    set(handles.UseImpactsListbox,'Visible','off')
    set(handles.CurrentFilePopup,'Visible','off')
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ProcessModePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProcessModePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PlotChannelsListbox.
function PlotChannelsListbox_Callback(hObject, eventdata, handles)
% hObject    handle to PlotChannelsListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PlotChannelsListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotChannelsListbox
n=get(hObject,'Value');
if isempty(n) |  n == 1
    n=handles.UseChannels;
    set(hObject,'Value',n);
else
    handles.UseChannels=n;
end
UpdateTimePlots(handles);
handles=UpdateFreqPlots(handles);
guidata(hObject,handles);

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
n=get(hObject,'Value');
if isempty(n)
    questdlg('You must use at least one impact!','OK');
else
    handles.UseImpacts=n;
end
set(handles.UseImpactsListbox,'Value',handles.UseImpacts);
UpdateTimePlots(handles);
handles=UpdateFreqPlots(handles);
guidata(hObject,handles)

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


% --- Executes on selection change in CurrentFilePopup.
function CurrentFilePopup_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentFilePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CurrentFilePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CurrentFilePopup
contents=cellstr(get(hObject,'String'));
FileName=contents{get(hObject,'Value')};
FileName=fullfile(handles.ImpactGui.InputDirectory,FileName);
handles.CurrentFile=FileName;
load(FileName,'-mat');
    % Check if impacts are negative, then "flip" positive
    handles.x=cell2mat(Data(1));
    [M,I]=max(abs(handles.x));
    if handles.x(I) < 0
        handles.x=-handles.x;
        Dir=dir2nbr(Header(1).Dir);
        Header.Dir=nbr2dir(-Dir);
    end
handles.x=cell2mat(Data(1));
handles.y=cell2mat(Data(2:end));
handles.Header=Header;
[handles.TrigIdx,handles.DIdx,handles.AIdx]=GuiImptrig(handles.x,handles.N,...
    handles.ImpactGui.TrigLevel,handles.ImpactGui.Pretrigger);
handles.UseImpacts=[1:length(handles.TrigIdx)];
set(handles.UseImpactsListbox,'String',{handles.UseImpacts'});
set(handles.UseImpactsListbox,'Value',handles.UseImpacts);
UpdateTimePlots(handles);
handles=UpdateFreqPlots(handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function CurrentFilePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CurrentFilePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
if isfield(handles,'h1') & ishghandle(handles.h1)
    close(handles.h1)
end
if isfield(handles,'h2') & ishghandle(handles.h2)
    close(handles.h2)
end
varargout{1} = '';
delete(hObject);

%=====================================================================
% User defined functions
%=====================================================================

function UpdateTimePlots(handles)
% Update the time plots of force and acceleration(s)
% Plot first layer
figure(handles.h1);
subplot(2,1,1)
t=makexaxis(handles.x,1/handles.fs);
plot(t,handles.x)
ylabel('Force','FontSize',8)
set(gca,'FontSize',8)
axis tight
hold on
% Put text with trigger number next to each impact. 
Ax=axis;
ylevel=0.9*Ax(4);              % text y level
for n = 1:length(handles.AIdx)
    text(t(handles.AIdx(n)),ylevel,int2str(n),'FontSize',7);
end
% Overlay with all triggered blocks in green
N=handles.ImpactGui.BlocksizeNumbers(handles.ImpactGui.BlocksizeNumber);
TrigIdx=handles.TrigIdx(handles.UseImpacts);
for n=1:length(TrigIdx)
    idx=TrigIdx(n):TrigIdx(n)+N-1;
    plot(t(idx),handles.x(idx),'g')
end
if handles.DIdx ~= 0
    for n=1:length(handles.DIdx)
        idx=handles.TrigIdx(handles.DIdx(n)):handles.TrigIdx(handles.DIdx(n))+N-1;
        plot(t(idx),handles.x(idx),'r')
    end
end
hold off
% Plot response(s)
subplot(2,1,2)
plot(t,handles.y(:,handles.UseChannels-1),'b')
ylabel('Response(s)','FontSize',8)
xlabel('Time [s]','FontSize',8)
axis tight
set(gca,'FontSize',8)
hold on
for n=1:length(TrigIdx)
    idx=TrigIdx(n):TrigIdx(n)+N-1;
    plot(t(idx),handles.y(idx,handles.UseChannels-1),'g')
end
if handles.DIdx ~= 0
    for n=1:length(handles.DIdx)
        idx=handles.TrigIdx(handles.DIdx(n)):handles.TrigIdx(handles.DIdx(n))+N-1;
        plot(t(idx),handles.y(idx,handles.UseChannels-1),'r')
    end
end
hold off

function FillUseFilesPopoup(handles)
% After an update in the input files settings, the popup menu with all
% imptime files in the current input directory with Prefix and StartNo,
% LastNo needs to be filled. If no files found, that is written to the
% popup menu
% Fill popup menu for imptime files in input directory
% Filelist=handles.ImpactGui.InputFileList);
% FileList={};
% for n = handles.ImpactGui.InputStartNo:handles.ImpactGui.InputLastNo
%     FileName=strcat(FilePrefix,int2str(n),'.imptime');
%     if exist(FileName,'file') == 2
%         FileList=[FileList; strcat(handles.ImpactGui.InputPrefix,int2str(n))];
%     end
% end
if length(handles.ImpactGui.InputFileList) == 1
    FileList={'No Files!'};
    set(handles.CurrentFilePopup,'Value',1);
end
set(handles.CurrentFilePopup,'String',handles.ImpactGui.InputFileList)



function handles = UpdateFreqPlots(handles)
% Update the Frequency domain plots

% Compute FRF(s) etc.
if ~isfield(handles,'flo')  % Manual methods
    [handles.H,handles.f, handles.C, handles.Tff] = Imp2FrfGui(handles.x,...
        handles.y,handles.fs,handles.N,handles.TrigIdx(handles.UseImpacts),...
        handles.DIdx,handles.ImpactGui.ForceWindowLength,...
        handles.ImpactGui.ExpWindowPercent,handles.ProcessType);
else    % optimization methods
    [handles.H,handles.f, handles.C, handles.Tff, handles.UseImpacts] = Imp2FrfGui(handles.x,...
        handles.y,handles.fs,handles.N,handles.TrigIdx(handles.UseImpacts),...
        handles.DIdx,handles.ImpactGui.ForceWindowLength,...
        handles.ImpactGui.ExpWindowPercent,handles.ProcessType,...
        handles.flo,handles.fhi);
end
figure(handles.h2);
% Plot FRF, Coherence, and Force Spectrum
subplot(4,1,[1 2])
semilogy(handles.f,abs(handles.H(:,handles.UseChannels-1)))
ylabel('Mag. FRF','FontSize',8)
axis tight
grid
set(gca,'FontSize',8)
subplot(4,1,3)
if isempty(handles.C)
    A=axis;
    plot(0,0)
    axis(A)
    text(0.1*(handles.f(end)-handles.f(1)),0.5,'Coherence undefined')
elseif handles.C == 0 
    A=axis;
    plot(0,0)
    axis(A)
    text(0.1*(handles.f(end)-handles.f(1)),0.5,'Coherence undefined')
else
    plot(handles.f,handles.C(:,handles.UseChannels-1))
    ylabel('Coherence','FontSize',8)
    axis tight
    grid
    ylim([0 1.05])
    set(gca,'FontSize',8)
end
subplot(4,1,4)
semilogy(handles.f,handles.Tff)
xlabel('Frequency [Hz]','FontSize',8)
ylabel('Force Spectrum','FontSize',8)
axis tight
grid
set(gca,'FontSize',8)



function handles = SaveFiles(handles)
% Save the FRFs, coherences, and force spectrum in handles.H, etc. 
% Also update handles fields for next file save
InHeader=handles.Header;
D=length(handles.y(1,:));
for k=1:D
    Header=InHeader(k+1);
    Header.RefDof=InHeader(1).Dof;
    Header.RefDir=InHeader(1).Dir;
    Header.Unit=strcat('(',InHeader(k+1).Unit,')/',InHeader(1).Unit);
    Header.FunctionType=4;
    Header.xIncrement=handles.f(2);
    Header.NoValues=length(handles.H);
    % Add special fields for impact FRFs
    Header.ExpWinEndPercent=handles.ImpactGui.ExpWindowPercent;
    if iscell(handles.UseImpacts)  % If ProcessMode is automatic
        Header.ImpactNumbers=handles.UseImpacts{k};
    else
        Header.ImpactNumbers=handles.UseImpacts;
    end
    Header.OrgFileName=handles.CurrentFile;
    Header.ProcessingMode=handles.ProcessType;
    FileName=strcat(handles.ImpactGui.OutputDirectory,'\',...
        handles.ImpactGui.OutputPrefix,'H',...
        int2str(handles.ImpactGui.OutputStartNo));
    Data=handles.H(:,k);
    if k == 1
        % First check if first file already exists, then ask if overwrite
        if exist(FileName,'file') == 2
            status=questdlg('File exists! Overwrite?');
        else
            status='YES';
        end
    end
    if strcmp(upper(status),'YES')
        if k == 1
            fprintf('Saving files to directory: %s\n',handles.ImpactGui.OutputDirectory);
        end
        fprintf('Saving file: %s\n',...
            strcat(handles.ImpactGui.OutputPrefix,'H',...
            int2str(handles.ImpactGui.OutputStartNo)));
        save(FileName,'Data','Header')
        Header.Unit='';
        Header.FunctionType=6;
        FileName=strcat(handles.ImpactGui.OutputDirectory,'\',...
            handles.ImpactGui.OutputPrefix,'C',...
            int2str(handles.ImpactGui.OutputStartNo));
        Data=handles.C(:,k);
        save(FileName,'Data','Header')
        Header.FunctionType=30;
        Header.Unit=strcat(InHeader(1).Unit,'s');
        Header.Dof=InHeader(1).Dof;
        Header.Dir=InHeader(1).Dir;
        Header=rmfield(Header,'RefDof');
        Header=rmfield(Header,'RefDir');
        FileName=strcat(handles.ImpactGui.OutputDirectory,'\',...
            handles.ImpactGui.OutputPrefix,'F',...
            int2str(handles.ImpactGui.OutputStartNo));
        Data=handles.Tff;
        save(FileName,'Data','Header')
        handles.ImpactGui.OutputStartNo=handles.ImpactGui.OutputStartNo+1;
    elseif strcmp(upper(status),'NO')
        ImpactContDlg('Title','Go back to Main GUI and change output directory')
    end
end
