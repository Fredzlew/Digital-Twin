function varargout = ImpactContDlg(varargin)
% IMPACTCONTDLG MATLAB code for ImpactContDlg.fig
%      IMPACTCONTDLG, by itself, creates a new IMPACTCONTDLG or raises the existing
%      singleton*.
%
%      H = IMPACTCONTDLG returns the handle to a new IMPACTCONTDLG or the handle to
%      the existing singleton*.
%
%      IMPACTCONTDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPACTCONTDLG.M with the given input arguments.
%
%      IMPACTCONTDLG('Property','Value',...) creates a new IMPACTCONTDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImpactContDlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImpactContDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2009-2014 by Anders Brandt
% Email: abra@iti.sdu.dk
% Version: 1.0 2014-07-15
% This file is part of ABRAVIBE Toolbox for NVA


% Last Modified by GUIDE v2.5 04-Jul-2014 15:12:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImpactContDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @ImpactContDlg_OutputFcn, ...
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


% --- Executes just before ImpactContDlg is made visible.
function ImpactContDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImpactContDlg (see VARARGIN)

% Choose default command line output for ImpactContDlg
handles.output = hObject;

set(handles.DlgText,'String',varargin{2})

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImpactContDlg wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImpactContDlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in ContinueButton.
function ContinueButton_Callback(hObject, eventdata, handles)
% hObject    handle to ContinueButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)
