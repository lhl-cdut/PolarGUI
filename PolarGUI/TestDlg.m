function varargout = TestDlg(varargin)
% TESTDLG MATLAB code for TestDlg.fig
%      TESTDLG, by itself, creates a new TESTDLG or raises the existing
%      singleton*.
%
%      H = TESTDLG returns the handle to a new TESTDLG or the handle to
%      the existing singleton*.
%
%      TESTDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTDLG.M with the given input arguments.
%
%      TESTDLG('Property','Value',...) creates a new TESTDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TestDlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TestDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TestDlg

% Last Modified by GUIDE v2.5 14-Dec-2018 12:28:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TestDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @TestDlg_OutputFcn, ...
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


% --- Executes just before TestDlg is made visible.
function TestDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TestDlg (see VARARGIN)

% Choose default command line output for TestDlg
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TestDlg wait for user response (see UIRESUME)
% uiwait(handles.TestDlg);


% --- Outputs from this function are returned to the command line.
function varargout = TestDlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
