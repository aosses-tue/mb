function varargout = learningGUI(varargin)
% LEARNINGGUI MATLAB code for learningGUI.fig
%      LEARNINGGUI, by itself, creates a new LEARNINGGUI or raises the existing
%      singleton*.
%
%      H = LEARNINGGUI returns the handle to a new LEARNINGGUI or the handle to
%      the existing singleton*.
%
%      LEARNINGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEARNINGGUI.M with the given input arguments.
%
%      LEARNINGGUI('Property','Value',...) creates a new LEARNINGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before learningGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to learningGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help learningGUI

% Last Modified by GUIDE v2.5 16-Jan-2015 01:03:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @learningGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @learningGUI_OutputFcn, ...
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


% --- Executes just before learningGUI is made visible.
function learningGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to learningGUI (see VARARGIN)

% Choose default command line output for learningGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes learningGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = learningGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnFactorial.
function btnFactorial_Callback(hObject, eventdata, handles)
% hObject    handle to btnFactorial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% https://www.youtube.com/watch?v=qAGiAq6HhYQ

n = str2num( get( handles.editNumber,'string' ) );

f = 1;
for i = 1:n
    f = f*i;
end

ff = num2str(f);

cc = get(handles.bDoAnalysis,'value');
set(handles.txtResult,'string',[ff ' - ' num2str(cc)]);

function editNumber_Callback(hObject, eventdata, handles)
% hObject    handle to editNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumber as text
%        str2double(get(hObject,'String')) returns contents of editNumber as a double


% --- Executes during object creation, after setting all properties.
function editNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bDoAnalysis.
function bDoAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to bDoAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bDoAnalysis


% --- Executes on button press in bDoSTFT.
function bDoSTFT_Callback(hObject, eventdata, handles)
% hObject    handle to bDoSTFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bDoSTFT
