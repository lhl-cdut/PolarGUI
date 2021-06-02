function varargout = FilterDlg(varargin)
% FILTERDLG MATLAB code for FilterDlg.fig
%      FILTERDLG, by itself, creates a new FILTERDLG or raises the existing
%      singleton*.
%
%      H = FILTERDLG returns the handle to a new FILTERDLG or the handle to
%      the existing singleton*.
%
%      FILTERDLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILTERDLG.M with the given input arguments.
%
%      FILTERDLG('Property','Value',...) creates a new FILTERDLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FilterDlg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FilterDlg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FilterDlg

% Last Modified by GUIDE v2.5 04-Jul-2018 13:21:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FilterDlg_OpeningFcn, ...
                   'gui_OutputFcn',  @FilterDlg_OutputFcn, ...
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


% --- Executes just before FilterDlg is made visible.
function FilterDlg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FilterDlg (see VARARGIN)

% Choose default command line output for FilterDlg
handles.output = hObject;


%--------------滤波器类型及算法选择控制
global g_ftype;
global g_falgr;

g_ftype = 1;
g_falgr = 1;

global isOpenFilter; %打开Filter标志

isOpenFilter = 0;

% set default value

global g_nsfre; %unit-ms

tfs = num2str(1000/g_nsfre); % Hz

set(handles.ed_fs,'string',tfs); 
set(handles.ed_lpfp,'string','200'); %通带截止频率
set(handles.ed_lpfs,'string','220'); %阻带起始频率


% band[fs1 fp1      fp2  fs2]%过度带=fp1-fs1
set(handles.ed_bpfs1,'string','0'); %阻带下限起始频率
set(handles.ed_bpfp1,'string',num2str(str2num(get(handles.ed_bpfs1,'string'))+20)); %通带下限截止频率
set(handles.ed_bpfp2,'string','200'); %通带上限截止频率
set(handles.ed_bpfs2,'string',num2str(str2num(get(handles.ed_bpfp2,'string'))+20)); %阻带上限起始频率


% set(ed_order,'string','20'); 

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FilterDlg wait for user response (see UIRESUME)
%uiwait(handles.bf_filter);


% --- Outputs from this function are returned to the command line.
function varargout = FilterDlg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%delete(handles.bf_filter);



% --- Executes on button press in rd_gauss.
function rd_gauss_Callback(hObject, eventdata, handles)
% hObject    handle to rd_gauss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rd_gauss


% --- Executes on button press in rd_kasier.
function rd_kasier_Callback(hObject, eventdata, handles)
% hObject    handle to rd_kasier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rd_kasier



% --- Executes on button press in rd_ham.
function rd_ham_Callback(hObject, eventdata, handles)
% hObject    handle to rd_ham (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rd_ham



function d_bpfl_Callback(hObject, eventdata, handles)
% hObject    handle to d_bpfl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_bpfl as text
%        str2double(get(hObject,'String')) returns contents of d_bpfl as a double



function d_bpfh_Callback(hObject, eventdata, handles)
% hObject    handle to d_bpfh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_bpfh as text
%        str2double(get(hObject,'String')) returns contents of d_bpfh as a double



function d_lpffre_Callback(hObject, eventdata, handles)
% hObject    handle to d_lpffre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_lpffre as text
%        str2double(get(hObject,'String')) returns contents of d_lpffre as a double



% --- Executes when selected object is changed in bf_filtertype.
function bf_filtertype_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bf_filtertype 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global g_ftype;
%strtype = get(hObject,'string');
strtype = get(eventdata.NewValue,'string');
switch strtype
    case 'LPFilter'
        g_ftype = 1;
    case 'BPFilter'
        g_ftype = 2;       
    otherwise
        msgbox('Please choose a filter type');
        g_ftype = 1;
end


% --- Executes when selected object is changed in bf_filteralg.
function bf_filteralg_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in bf_filteralg 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global g_falgr;
%stralgr = get(hObject,'string');
stralgr = get(eventdata.NewValue,'string');
switch stralgr 
    
    case 'Hamming'
        g_falgr = 1;
        
    case 'Gauss'
        g_falgr = 2;
       
    case 'Kasier'
        g_falgr = 3;
        
    case 'Hanning'
        g_falgr = 4;
        
    case 'Blackman'
        g_falgr = 5;
        
    case 'Butterworth'
        g_falgr = 6;

    otherwise
        msgbox('please choose a filter algrithm');
        g_falgr = 1;
end

% --- Executes on button press in b_filterok.
function b_filterok_Callback(hObject, eventdata, handles)
% hObject    handle to b_filterok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%读取地震数据
global g_ntrace;
global g_nslen;
global g_nsfre;
% global g_oriseisdata


global g_filterdata;


%--------------滤波器类型及算法选择控制
global g_ftype;
global g_falgr;

global isOpenFilter; %打开Filter标志

isOpenFilter = 1;
global isButter;
isButter = 0;
switch g_ftype
    case 1  %FIR低通滤波
        g_ftype = 1;
        fs = str2num(get(handles.ed_fs,'string'));% transfer the sample rate in Hz
        fcl = str2num(get(handles.ed_lpfp,'string')); %the edge frequency of passband in Hz
        fch = str2num(get(handles.ed_lpfs,'string')); %the edge frequency of stopband in Hz
        wp = 2*pi*fcl/fs; %normalization angular frequency for the passband
        ws = 2*pi*fch/fs; %normalization angular frequency for the stopband
        wn = (wp+ws)/(2*pi); % cut off frequency;
        Bt = ws-wp; % transition band
        N0 = ceil(6.2*pi/Bt)+1; %filter orders
        N = N0+mod(N0+1,2);
        if N>51 
            N=51; % Keep order <= 51
        end 
        
        winHan = hanning(N); %returns the N-point symmetric Hanning window
        winHam = hamming(N);%使用hamming窗函数
        winBlkman = blackman(N);%使用blackman窗函数
        winGauss = gausswin(N);
        [n,Wn,beta,ftype] = kaiserord([20 25],[1 0],[0.01 0.01],100);
        winKs = kaiser(n+1,beta);%使用kaiser窗函数
        
        switch g_falgr
            case 1 %Haming
                g_falgr = 1;
                [b,a] = fir1(N-1,wn,'low',winHam);
                
            case 2 %Gauss
                g_falgr = 1;
                [b,a] = fir1(N-1,wn,'low',winGauss);
                
            case 3 %Kasier
                g_falgr = 1;
                [b,a] = fir1(n,Wn/pi,'low',winKs,'noscale');
                
            case 4 % Hanning
                g_falgr = 1;
                [b,a] = fir1(N-1,wn,'low',winHan);
                
            case 5 %Blackman
                g_falgr = 1;
                [b,a] = fir1(N-1,wn,'low',winBlkman);
                
            case 6 %Butterworth
                g_falgr = 1;
                isButter = 1;
                
            otherwise
                msgbox('Not choose filter algrithm');
        end
        if (isButter == 1) %only for earthquake data sac or mseed
            isButter = 0;
            dt = g_nsfre/1000; %ms-s
            for k = 1:g_ntrace %所有道滤波
%                 g_filterdata(:,k) = buttern_low(g_seisdata(:,k),N,wn,dt);

            trow = buttern_low(g_filterdata(:,k)',6,(fcl+fch)/2,dt);
            g_filterdata(:,k) = trow';
            %default fch = 1Hz only for earthquake data sac or mseed
            
            end
        else

            % g_filterdata = filter(b,a,1,g_filterdata);
             g_filterdata = filtfilt(b,a,g_filterdata);

            %      g_filterdata = Kaiser(g_seisdata,0,1,uint8(1/g_nsfre));
        end

    case 2  %带通滤波
        g_ftype = 1;
        Fs = str2num(get(handles.ed_fs,'string'));
        
        fp1 = str2num(get(handles.ed_bpfp1,'string')); %通带下限截止频率
        fp2 = str2num(get(handles.ed_bpfp2,'string')); %通带上限截止频率
        fs1 = str2num(get(handles.ed_bpfs1,'string')); %阻带下限截止频率
        fs2 = str2num(get(handles.ed_bpfs2,'string')); %阻带上限截止频率
        
        wp1 = 2*pi*fp1/Fs;%将通带下限截止频率转换为数字滤波器频率
        wp2 = 2*pi*fp2/Fs;%将通带上限截止频率转换为数字滤波器频率
        ws1 = 2*pi*fs1/Fs;%将通带下限截止频率转换为数字滤波器频率
        ws2 = 2*pi*fs2/Fs;%将通带上限截止频率转换为数字滤波器频率
        
        Bt = wp1-ws1; %过渡带
        N0 = ceil(6.2*pi/Bt)+1;
        N = N0+mod(N0+1,2);
        if N>101
            N=101;  % Keep order < 101
        end
        wn = [(wp1+ws1)/(2*pi),(wp2+ws2)/(2*pi)];
        winHan = hanning(N);%使用hanning窗函数
        winHam = hamming(N);%使用hamming窗函数
        winBlkman = blackman(N);%使用blackman窗函数
        winGauss = gausswin(N);
        %设过渡带宽度为5Hz
        [n,Wn,beta,ftype] = kaiserord([10 15 20 25],[0 1 0],[0.01 0.01 0.01],100);%求阶数n以及参数beta
        winKs = kaiser(n+1,beta);%使用kaiser窗函数
    
        
        switch g_falgr
            case 1 %Haming
                g_falgr = 1;
                [b,a] = fir1(N-1,wn,'bandpass',winHam);
                
            case 2 %Gauss
                g_falgr = 1;
                [b,a] = fir1(N-1,wn,'bandpass',winGauss);
                
            case 3 %Kasier
                g_falgr = 1;
                [b,a] = fir1(n,Wn,winKs,'bandpass','noscale');
                
            case 4 % Hanning
                g_falgr = 1;
                [b,a] = fir1(N-1,wn,'bandpass',winHan);
                
            case 5 %Blackman
                g_falgr = 1;
                [b,a] = fir1(N-1,wn,'bandpass',winBlkman);
                
            case 6 %Butterworth
                g_falgr = 1;
                isButter = 1;
    
            otherwise
                msgbox('Not choose filter algrithm');
        end
        
        if (isButter == 1) %only for earthquake data sac or mseed file
            isButter = 0;
            dt = g_nsfre/1000; %ms-s
            for k = 1:g_ntrace
               trow = buttern_filter(g_filterdata(:,k)',6,(fs1+fp1)/2,(fs2+fp2)/2,dt);
               g_filterdata(:,k) = trow';
            %Default fp1=1/50,fp2=1/10; %only for earthquake data sac or mseed file
            end
          
        else
%             g_filterdata = filter(b,a,1,g_filterdata);            
            g_filterdata = filtfilt(b,a,g_filterdata);
        end
        
    otherwise
        msgbox('Not choose filter tpye');
        g_falgr = 1;

end

%isOpenFilter = 1;
%hold off;
% hmshandles=guihandles(MSLocation);
% uiresume(hmshandles.MSLocation);%用户操作完成  
%close(gcf);%关闭子窗口
delete(handles.bf_filter);


% --- Executes on button press in b_filterexit.
function b_filterexit_Callback(hObject, eventdata, handles)
% hObject    handle to b_filterexit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global isOpenFilter; %打开Filter标志
isOpenFilter = 0; %打开Filter等待处理标志
delete(handles.bf_filter);




function ed_lpfp_Callback(hObject, eventdata, handles)
% hObject    handle to ed_lpfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_lpfp as text
%        str2double(get(hObject,'String')) returns contents of ed_lpfp as a double


% --- Executes during object creation, after setting all properties.
function ed_lpfp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_lpfp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_lpfs_Callback(hObject, eventdata, handles)
% hObject    handle to tx_lpfs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tx_lpfs as text
%        str2double(get(hObject,'String')) returns contents of tx_lpfs as a double


% --- Executes during object creation, after setting all properties.
function tx_lpfs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_lpfs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rd_hann.
function rd_hann_Callback(hObject, eventdata, handles)
% hObject    handle to rd_hann (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rd_hann



function ed_fs_Callback(hObject, eventdata, handles)
% hObject    handle to ed_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_fs as text
%        str2double(get(hObject,'String')) returns contents of ed_fs as a double


% --- Executes during object creation, after setting all properties.
function ed_fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ed_lpfs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_lpfs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_bpfp1_Callback(hObject, eventdata, handles)
% hObject    handle to ed_bpfp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_bpfp1 as text
%        str2double(get(hObject,'String')) returns contents of ed_bpfp1 as a double


% --- Executes during object creation, after setting all properties.
function ed_bpfp1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_bpfp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_bpfp2_Callback(hObject, eventdata, handles)
% hObject    handle to ed_bpfp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_bpfp2 as text
%        str2double(get(hObject,'String')) returns contents of ed_bpfp2 as a double


% --- Executes during object creation, after setting all properties.
function ed_bpfp2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_bpfp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_bpfs1_Callback(hObject, eventdata, handles)
% hObject    handle to ed_bpfs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_bpfs1 as text
%        str2double(get(hObject,'String')) returns contents of ed_bpfs1 as a double


% --- Executes during object creation, after setting all properties.
function ed_bpfs1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_bpfs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ed_bpfs2_Callback(hObject, eventdata, handles)
% hObject    handle to ed_bpfs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_bpfs2 as text
%        str2double(get(hObject,'String')) returns contents of ed_bpfs2 as a double


% --- Executes during object creation, after setting all properties.
function ed_bpfs2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_bpfs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
