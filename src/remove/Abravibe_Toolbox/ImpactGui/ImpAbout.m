function varargout = ImpAbout(varargin)
% IMPABOUT MATLAB code for ImpAbout.fig
%      IMPABOUT, by itself, creates a new IMPABOUT or raises the existing
%      singleton*.
%
%      H = IMPABOUT returns the handle to a new IMPABOUT or the handle to
%      the existing singleton*.
%
%      IMPABOUT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPABOUT.M with the given input arguments.
%
%      IMPABOUT('Property','Value',...) creates a new IMPABOUT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImpAbout_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImpAbout_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (c) 2009-2014 by Anders Brandt
% Email: abra@iti.sdu.dk
% Version: 1.0 2014-07-15
% This file is part of ABRAVIBE Toolbox for NVA


% Last Modified by GUIDE v2.5 11-Jul-2014 16:04:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImpAbout_OpeningFcn, ...
                   'gui_OutputFcn',  @ImpAbout_OutputFcn, ...
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


% --- Executes just before ImpAbout is made visible.
function ImpAbout_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImpAbout (see VARARGIN)

% Choose default command line output for ImpAbout
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImpAbout wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImpAbout_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
