%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% Polarization_GUI.m
% Polarization analysis for 3C seismic recordings adopting 4 algorithms.
%
% PROGRAMMER:
% Huailiang Li
%
% Last revision date:
% 09 Dec 2019
% Last modified by ()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Polarization_GUI

close(gcf);
close all;
clear;
clc;
% Include all sub-paths
addpath(genpath(pwd));

global g_backColor;
% g_backColor = [0.11 0.47 0.53];
g_backColor = [0.22 0.34 0.42];

% Curve color
gcolorR1 = [0.5 0 0];
gcolorR = [0.72 0.15 0.2];
gcolorG = [0.25 0.45 0.2];
gcolorB = [0 0.44 0.74];


global g_ntrace;
global g_nslen;
global g_nsfre;
global g_seisdata;
global g_filterdata;
global g_format;
global g_hf;

global g_Pcursorb;   % Begin of p-wave cut window for cursor handle
global g_Pcursore;   % End of p-wave cut window for cursor handle
global g_Scursorb;   % Begin of p-wave cut window for cursor handle
global g_Scursore;   % End of p-wave cut window for cursor handle

global g_ptb;
global g_pte;
g_ptb  = zeros(1,1);   % Begin of p-wave cut window
g_pte  = zeros(1,1);   % End of p-wave cut window
global g_stb;
global g_ste;
g_stb  = zeros(1,1);   % Begin of p-wave cut window
g_ste  = zeros(1,1);   % End of p-wave cut window


g_ntrace = double(zeros(0));
g_nslen = double(zeros(0));
g_nsfre = double(zeros(0));

global isOpenFilter;
isOpenFilter = 0;

global isOpenSAC;
isOpenSAC = 0;
global g_polarimod;  % Polarization method
g_polarimod = 1;

global bg_polarizmehtod;

global isOpennum;
isOpennum = 0;

global g_PrinFreP;
g_PrinFreP = 100;
global g_PrinFreS;
g_PrinFreS = 100;

global IsOpenR3D;
IsOpenR3D = 0;


global hfgui;

%===================================================================================
% Main Window Initial
%===================================================================================
hfgui = figure(...
    'Name', 'PolarizationGUI',...
    'NumberTitle', 'off',...
    'toolbar','none',...
    'Units', 'Normalized', ...
    'Menu', 'none',...
    'Color', g_backColor,...  %Color
    'Position', [0.02 0.01 0.92 0.99]);

% figuresize = get(hfgui,'Position');
% screensize = get(0,'screensize');

% Openfile Button
hbtn_openfile = uicontrol(...
    hfgui,...
    'Style', 'pushbutton',...
    'Units', 'Normalized',...
    'Position', [0.00 0.92 0.05 0.08],...
    'Callback',@(x,y)pb_openfile_Callback);
imgopenfile = imread('openfile.jpg');
set(hbtn_openfile,'CDATA',imgopenfile);
set(hbtn_openfile,'foregroundcolor',g_backColor,'backgroundcolor',g_backColor);

% Exit Button
hbtn_exit = uicontrol(...
    hfgui,...
    'Style', 'pushbutton',...
    'Units', 'Normalized',...
    'Position', [0.05 0.92 0.05 0.08],...
    'Callback',@(~,~)pb_exit_Callback);
imgexit = imread('exit.png');
set(hbtn_exit,'CDATA',imgexit);
set(hbtn_exit,'foregroundcolor',g_backColor,'backgroundcolor',g_backColor);

%--------------------------------------------------------------------------
% Polarization Method Choose
bg_polarizmehtod = uibuttongroup(hfgui,...
    'unit','Normalized',...
    'position',[0.00 0.78 0.10 0.14],...
    'title','Polarization Method',...
    'TitlePosition','centertop',...
    'fontsize',12,...
    'ShadowColor',g_backColor,...
    'backgroundcolor',g_backColor,...
    'foregroundcolor','w');

rd_covar1 = uicontrol(bg_polarizmehtod,...
    'style','radiobutton',...
    'unit','Normalized',...
    'string','EVD-CM1',...
    'Tag','EVD-CM1',...
    'Value',0,...
    'fontsize',12,...
    'HandleVisibility','off',...
    'position',[0.10 0.75 1.00 0.25],...
    'backgroundcolor',g_backColor,...
    'foregroundcolor','w');

rd_covar2 = uicontrol(bg_polarizmehtod,...
    'style','radiobutton',...
    'unit','Normalized',...
    'position',[0.10 0.50 1.00 0.25],...
    'string','EVD-CM2',...
    'Tag','EVD-CM2',...
    'fontsize',12,...
    'HandleVisibility','off',...
    'backgroundcolor',g_backColor,...
    'foregroundcolor','w');

rd_SVD = uicontrol(bg_polarizmehtod,...
    'style','radiobutton',...
    'unit','Normalized',...
    'position',[0.10 0.25 1.00 0.25],...
    'string','PCA-SVD',...
    'Tag','PCA-SVD',...
    'fontsize',12,...
    'HandleVisibility','off',...
    'backgroundcolor',g_backColor,...
    'foregroundcolor','w');

rd_AnaSig = uicontrol(bg_polarizmehtod,...
    'style','radiobutton',...
    'unit','Normalized',...
    'position',[0.1 0.00 1.00 0.25],...
    'string','EVD-ASM',...
    'Tag','EVD-ASM',...
    'fontsize',12,...
    'HandleVisibility','off',...
    'backgroundcolor',g_backColor,...
    'foregroundcolor','w');

set(bg_polarizmehtod,'Visible','on');
set(bg_polarizmehtod,'SelectionChangeFcn',@btp_SelFcn);
%--------------------------------------------------------------------------
% Filter Button
hbtn_filter = uicontrol(...
    'Style', 'pushbutton',...
    'Units', 'Normalized',...
    'Position', [0.00 0.75 0.05 0.03],...
    'string','Filter',...
    'fontsize',16,...
    'Callback',@(~,~)pb_filter_Callback);
set(hbtn_filter,'foregroundcolor','w','backgroundcolor',g_backColor);

% Rotate3D Button
hbtn_r3d = uicontrol(...
    'Style', 'pushbutton',...
    'Units', 'Normalized',...
    'Position', [0.05 0.75 0.05 0.03],...
    'string','R3DOff',...
    'fontsize',16,...
    'Callback',@(~,~)pb_r3d_Callback);
set(hbtn_r3d,'foregroundcolor','w','backgroundcolor',g_backColor);
tipstr = sprintf('Rotate 3D Off or On \n Please close the Rotate3D mode to select the waveform using cursor');
set(hbtn_r3d,'Tooltipstring',tipstr);

% Frequency edit/text
htxt_p = uicontrol(...
    'Style', 'text',...
    'Units', 'Normalized',...
    'Position', [0.00 0.72 0.01 0.03],...
    'string','P',...
    'fontsize',12);
set(htxt_p,'foregroundcolor','w','backgroundcolor',g_backColor);

hedit_p = uicontrol(...
    'Style', 'edit',...
    'Units', 'Normalized',...
    'Position', [0.01 0.72 0.03 0.03]);
set(hedit_p,'foregroundcolor','w','backgroundcolor',g_backColor);

htxt_s = uicontrol(...
    'Style', 'text',...
    'Units', 'Normalized',...
    'Position', [0.04 0.72 0.01 0.03],...
    'string','S',...
    'fontsize',12);
set(htxt_s,'foregroundcolor','w','backgroundcolor',g_backColor);

hedit_s = uicontrol(...
    'Style', 'edit',...
    'Units', 'Normalized',...
    'Position', [0.05 0.72 0.03 0.03]);
set(hedit_s,'foregroundcolor','w','backgroundcolor',g_backColor);

htxt_hz = uicontrol(...
    'Style', 'text',...
    'Units', 'Normalized',...
    'Position', [0.08 0.72 0.02 0.03],...
    'string','Hz',...
    'fontsize',12);
set(htxt_hz,'foregroundcolor','w','backgroundcolor',g_backColor);

% Default Principal Frequency
set(hedit_p,'string',num2str(g_PrinFreP));
set(hedit_s,'string',num2str(g_PrinFreS));


axes_swave = axes('Position', [0.14 0.90 0.85 0.09],...
    'Units','Normalized',...
    'Box','on');
set(axes_swave,'xTick',[],'ytick',[]);
set(axes_swave,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor); 

axes_cutpwave = axes('Position', [0.14 0.78 0.40 0.09],...
    'Units','Normalized',...
    'Box','on');
set(axes_cutpwave,'xTick',[],'ytick',[]);
set(axes_cutpwave,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

axes_cutswave = axes('Position', [0.59 0.78 0.40 0.09],...
    'Units','Normalized',...
    'Box','on');
set(axes_cutswave,'xTick',[],'ytick',[]);
set(axes_cutswave,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);


axes_P3dtraj = axes('Position', [0.03 0.38 0.27 0.32],...
    'Units','Normalized',...
    'Box','on');
set(axes_P3dtraj,'xTick',[],'ytick',[]);
set(axes_P3dtraj,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);
 
axes_S3dtraj = axes('Position', [0.03 0.02 0.27 0.32],...
    'Units','Normalized',...
    'Box','on');
set(axes_S3dtraj,'xTick',[],'ytick',[]);
set(axes_S3dtraj,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);


% ROW 1
axes_Pazirose = axes('Position', [0.34 0.53 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_Pazirose,'xtick',[],'ytick',[],'xcolor',g_backColor,'ycolor',g_backColor);
set(axes_Pazirose,'color',g_backColor);

axes_Pincirose = axes('Position', [0.51 0.53 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_Pincirose,'xtick',[],'ytick',[],'xcolor',g_backColor,'ycolor',g_backColor);
set(axes_Pincirose,'color',g_backColor);

axes_Sazirose = axes('Position', [0.68 0.53 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_Sazirose,'xtick',[],'ytick',[],'xcolor',g_backColor,'ycolor',g_backColor);
set(axes_Sazirose,'color',g_backColor);

axes_Sincirose = axes('Position', [0.85 0.53 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_Sincirose,'xtick',[],'ytick',[],'xcolor',g_backColor,'ycolor',g_backColor);
set(axes_Sincirose,'color',g_backColor);

% ROW 2
axes_Pazimuth = axes('Position', [0.34 0.30 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_Pazimuth,'xTick',[],'ytick',[]);
set(axes_Pazimuth,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

axes_Pincidence = axes('Position', [0.51 0.30 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_Pincidence,'xTick',[],'ytick',[]);
set(axes_Pincidence,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

axes_Sazimuth = axes('Position', [0.68 0.30 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_Sazimuth,'xTick',[],'ytick',[]);
set(axes_Sazimuth,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

axes_Sincidence = axes('Position', [0.85 0.30 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_Sincidence,'xTick',[],'ytick',[]);
set(axes_Sincidence,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

% ROW 3
axes_PDlp = axes('Position', [0.34 0.04 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_PDlp,'xTick',[],'ytick',[]);
set(axes_PDlp,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

axes_PDpp = axes('Position', [0.51 0.04 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_PDpp,'xTick',[],'ytick',[]);
set(axes_PDpp,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

axes_SDlp = axes('Position', [0.68 0.04 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_SDlp,'xTick',[],'ytick',[]);
set(axes_SDlp,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

axes_SDpp = axes('Position', [0.85 0.04 0.12 0.20],...
    'Units','Normalized',...
    'Box','on');
set(axes_SDpp,'xTick',[],'ytick',[]);
set(axes_SDpp,'color',g_backColor,'XColor',g_backColor,'YColor',g_backColor);

rotate3d Off;
%===================================================================================
% Read seismic data
%===================================================================================

%--- Executes on button press in pb_openfile.
    function pb_openfile_Callback()
        % hObject    handle to pb_openfile (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        % g_handles = handles;
        
        vseisdata = double(zeros(0));
        
        if(isOpennum == 1)
            cla;
            cla reset;
            isOpennum = 0;
            
        end
        
        % Data order->X/Y/Z
        [filename, pathname] = uigetfile({'*.*';'*.mat';'*.asc';'*.csv';'*.dat';...
            '*.sg2';'*.sgy';'*.sac';'*.xls';'*.xlsx'}, 'load seismic data','MultiSelect','on'); % Choose a file
        % 'MultiSelect','on'-----only for opening SAC file
        if isequal(filename,0)   % Is choosing a file?
            msgbox('Choose no file');
            return;
        else
            pathfile=fullfile(pathname, filename);  % Obtaining file path
            if (iscell(filename)) % for Multi sac files reading
                [~,~,exc]=fileparts(pathfile{1,1});
                filenums = size(filename,2);  % for Multi files reading
            else
                [~,~,exc]=fileparts(pathfile);% for Multi files reading
            end
        end
        
        
        switch exc
            case '.dat'
                
                %         msgbox('pleas waiting');
                isOpennum = 1;
                g_format = '.dat';
                
            case {'.sg2','.SG2'}
                %--------------------------------------------------------------------------
                % using read_sg2 to load seismic data
                %function [ntraces, nsamplenums, nsamplefre, sg2data] = read_sg2(pathfile)
                %--------------------------------------------------------------------------
                %         [filename, pathname] = uigetfile('*.*', 'load seismic data');
                %         if isequal(filename,0)
                %             msgbox('Choose no ');
                %         else
                %             pathfile=fullfile(pathname, filename);
                %         end
                
                %         [vtrace,vslen,vsfre,vseisdata] = read_sg2(pathfile);
                
                [filehdr,trchdr,trctxt,data] = seg2read(pathfile, 'want','filehdr,trchdr,trctxt,data');
                vslen = trchdr(1).number_of_samples;
                vtrace = filehdr.number_of_traces_in_file;
                vsfre = str2num(trctxt(1).sample_interval)*1000; %s->ms
                vseisdata = data;
                isOpennum = 1;
                g_format = '.sg2';
                %--------------------------------------------------------------------------
                
            case {'.sgy','.SGY'}
                
                if(iscell(filename))
                    vtrace = 0;
                    for f = 1:1:filenums   % For separate X/Y/Z Files, Must be same length
                        [sdata,vslen,vsfre,strace,Theader] = altreadsegy(pathfile{1,f},'textheader','yes','fpformat','ieee');
                        
                        for m = 1:1:strace
                            vseisdata(:,((m-1)*3+f)) = sdata(:,m);
                        end
                        vtrace = vtrace + strace;
                    end
                    if (vtrace > 24)  % Cut 24 traces
                        vseisdata = vseisdata(:,(1:24));
                        vtrace = 24;
                    end
                else
                    [vseisdata,vslen,vsfre,vtrace,Theader] = altreadsegy(pathfile,'textheader','yes','fpformat','ieee');
                end
                
                isOpennum = 1;
                g_format = '.sgy';
                
                
            case '.mat'
                
                sydata = importdata(pathfile);
                vseisdata(:,1) = sydata.X;
                vseisdata(:,2) = sydata.Y;
                vseisdata(:,3) = sydata.Z;
                vtrace = 3;
                vslen = length(vseisdata(:,3));
                vsfre = 1/6;
                isOpennum = 1;
                g_format = '.mat';
                
            case '.asc'
                
                [X,Y,Z] = textread(pathfile,'%f,%f,%f','headerlines',1);
                vseisdata(:,1) = X;
                vseisdata(:,2) = Y;
                vseisdata(:,3) = Z;
                vtrace = 3;
                vslen = length(vseisdata(:,3));
                vsfre = 1/6 ;
                isOpennum = 1;
                g_format = '.asc';
                
            case '.csv'
                
                numerric = csvread(pathfile,1);
                [m,n]=size(numerric);
                %vseisdata = numerric(:,11:end);
                vseisdata = numerric;
                vtrace = n;
                vslen = m;
                vsfre = 1 ;
                isOpennum = 1;
                g_format = '.csv';
                
            case {'.sac','.SAC'}
                
                isOpenSAC = 1;
                if(iscell(filename))
                    
                    for f = 1:1:filenums
                        [vseisdata(:,f),tlen,headsac] = rdsac(pathfile{1,f});
                    end
                else
                    [vseisdata,tlen,headsac] = rdsac(pathfile);
                end
                
                vtrace = size(vseisdata,2);
                vslen = size(vseisdata,1);
                vsfre = headsac.DELTA*1000 ; %0.005s->ms
                isOpennum = 1;
                g_format = '.sac';
                
            case {'.mseed','.MSEED'}
                %         [X,I] = rdmseed(pathfile);
                %         X = rdmseed(pathfile);
                
                % test if the file is multiplexed or a single channel
                %	- to get the list of channels in a multiplexed file:
                
                
                if(iscell(filename))
                    
                    for f = 1:1:filenums
                        [X,I] = rdmseed(pathfile{1,f});
                        un = unique(cellstr(char(X.ChannelFullName)));
                        nc = numel(un); % traces nums
                        for i = 1:nc
                            k = I(i).XBlockIndex;
                            d(:,i+(f-1)*3) = cat(1,X(k).d);
                        end
                    end
                    
                    vtrace = nc*f;
                else
                    [X,I] = rdmseed(pathfile);
                    un = unique(cellstr(char(X.ChannelFullName)));
                    nc = numel(un); % traces nums
                    for i = 1:nc
                        k = I(i).XBlockIndex;
                        d(:,i) = cat(1,X(k).d);
                    end
                    vtrace = nc;
                end
                vseisdata = d;
                vslen = size(d,1);
                vsfre = 1/unique(cat(1,X.SampleRate))*1000; %Hz->s->ms
                isOpennum = 1;
                g_format = '.mseed';
                
            case {'.xls','.xlsx'} %for yangning data
                
                
                if(iscell(filename))
                    
                    for f = 1:1:filenums
                        fulldata = (xlsread(pathfile{1,f}))';
                        tnums = size(fulldata,2);% traces of signle file
                        for m = 1:1:tnums
                            vseisdata(:,tnums*(f-1)+m) = fulldata(:,m); %track one trace
                        end
                    end
                else
                    vseisdata = (xlsread(pathfile))';
                end
                
                vtrace = size(vseisdata,2);
                vslen = size(vseisdata,1);
                vsfre = 2; %2ms  for yang forward data
                isOpennum = 1;
                g_format = '.xls';
                
                
            case {'.txt','.TXT'} %for mengxiaobo-Fracturing log data
                
                
                if(iscell(filename))
                    
                    for f = 1:1:filenums
                        [tdata,vdata] = textread(pathfile{1,f},'%f %f');
                        vseisdata(:,f) = vdata; %track one trace
                    end
                else
                    [~,vseisdata] = textread(pathfile,'%f %f');
                end
                
                vtrace = size(vseisdata,2);
                vslen = size(vseisdata,1);
                vsfre = 0.5; %0.5ms/2kHz for mengxiaobo fracturing log data
                isOpennum = 1;
                g_format = '.txt';
                
            case {'.lvm'}
                d = lvm_import(pathfile,2);
                %                 vseisdata(:,1) = d.Segment1.data(:,4);
                %                 vseisdata(:,2) = d.Segment1.data(:,10);
                %                 vseisdata(:,3) = d.Segment1.data(:,5);
                vseisdata(:,1) = d.Segment1.data(:,7);
                vseisdata(:,2) = d.Segment1.data(:,4);
                vseisdata(:,3) = d.Segment1.data(:,6);
                vtrace = 8;
                vslen = size(vseisdata,1);
                vsfre = 1/6; %0.005s->ms
                isOpennum = 1;
                g_format = '.lvm';
                
            otherwise
                
                msgbox('Please choose a right format file');
                return;
                
        end
        
        g_ntrace = vtrace;
        g_nslen = vslen;
        g_nsfre = vsfre;
        
        %-------------------------------------
        % for making Figure in the paper, Please delete this part when using
%         g_seisdata = vseisdata(2500:3999,:); % for synthetic data-yang
%         g_seisdata = vseisdata(343000:345999,:); % for mengxiaobo fracturing log data
%         g_seisdata = vseisdata(5000:end,:); % for earthquake data with sac
%         g_seisdata = vseisdata(36800:38000,:);  % for active source microseismic data with lvm
        g_seisdata = vseisdata;
        g_filterdata = g_seisdata;
        g_nslen = length(g_seisdata);
        %-------------------------------------
        
        yff = fft(g_filterdata(:,3),g_nslen); % dataZ
        mag = abs(yff);
        fs = round(1000/g_nsfre); %g_nsfre (ms)->Hz
        [~,ind] = max(mag);
        g_PrinFreP = ind*fs/g_nslen;
        g_PrinFreS = g_PrinFreP;
        set(hedit_p,'string',num2str(g_PrinFreP));
        set(hedit_s,'string',num2str(g_PrinFreS));
        
        %--------------------------------------------------------------------------
        ha = axes_swave;
        Drawwave(ha,g_filterdata);
        
        g_hf = hfgui;
        
        if (isOpennum == 1)
            
            %-------------------------------------------------
            % For natural earthquake data (sac)
            %     g_Pcursorb = CreateCursor(g_hf,1,ha);
            %     SetCursorLocation(g_hf,g_Pcursorb, 11150+1000);
            %     g_Pcursore = CreateCursor(g_hf,2,ha);
            %     SetCursorLocation(g_hf,g_Pcursore, 11150+1600);
            %
            %     g_Scursorb = CreateCursor(g_hf,3,ha);
            %     SetCursorLocation(g_hf,g_Scursorb, 11150+3310);
            %     g_Scursore = CreateCursor(g_hf,4,ha);
            %     SetCursorLocation(g_hf,g_Scursore, 11150+3750);
            %-------------------------------------------------
            
            %-------------------------------------------------
            % For sythetic data-yang
            g_Pcursorb = CreateCursor(g_hf,1,ha);
            SetCursorLocation(g_hf,g_Pcursorb, g_nslen/4);
            g_Pcursore = CreateCursor(g_hf,2,ha);
            SetCursorLocation(g_hf,g_Pcursore, g_nslen/3);
            
            g_Scursorb = CreateCursor(g_hf,3,ha);
            SetCursorLocation(g_hf,g_Scursorb, g_nslen/2);
            g_Scursore = CreateCursor(g_hf,4,ha);
            SetCursorLocation(g_hf,g_Scursore, g_nslen); 
        end
        
        g_ptb = GetCursorLocation(g_hf,g_Pcursorb);
        g_pte = GetCursorLocation(g_hf,g_Pcursore);
        
        g_stb = GetCursorLocation(g_hf,g_Scursorb);
        g_ste = GetCursorLocation(g_hf,g_Scursore);
        
        Pdata = g_filterdata(g_ptb:g_pte,:);
        Sdata = g_filterdata(g_stb:g_ste,:);
        
        DrawCutwave(axes_cutpwave,Pdata);
        DrawCutwave(axes_cutswave,Sdata);
        
        delt = g_nsfre/1000; % ms->s
        ttotP = (g_pte-g_ptb+1)*delt;
        ttotS = (g_ste-g_stb+1)*delt;
        
        DrawAllGraph();
    end
%--------------------------------------------------------------------------

% --- Executes on button press in pb_exit.
    function pb_exit_Callback()
        % hObject    handle to pb_exit (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        
        close(gcf);
        close all;
        clear;
        clc;
    end
% --- Executes on button press in pb_filter.
    function pb_filter_Callback()
        % hObject    handle to pb_filter (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)       
        
        if (isOpennum == 1)
            
            g_ptb = GetCursorLocation(g_hf,g_Pcursorb);
            g_pte = GetCursorLocation(g_hf,g_Pcursore);
            g_stb = GetCursorLocation(g_hf,g_Scursorb);
            g_ste = GetCursorLocation(g_hf,g_Scursore);
            
            hflterhandles = guihandles(FilterDlg);
            uiwait(hflterhandles.bf_filter);  
            
            Pdata = g_filterdata(g_ptb:g_pte,:);
            Sdata = g_filterdata(g_stb:g_ste,:);
            
            if(isOpenFilter == 1)
                
                isOpenFilter = 0;
                
                DeleteCursor(g_Pcursorb);
                DeleteCursor(g_Pcursore);
                DeleteCursor(g_Scursorb);
                DeleteCursor(g_Scursore);
                
                ha = axes_swave;
                cla(ha);
                cla(axes_cutpwave);
                cla(axes_cutswave);
                
                Drawwave(ha,g_filterdata);
                
                g_Pcursorb = CreateCursor(g_hf,1,ha);
                SetCursorLocation(g_hf,g_Pcursorb, g_ptb);
                g_Pcursore = CreateCursor(g_hf,2,ha);
                SetCursorLocation(g_hf,g_Pcursore, g_pte);
                g_Scursorb = CreateCursor(g_hf,3,ha);
                SetCursorLocation(g_hf,g_Scursorb, g_stb);
                g_Scursore = CreateCursor(g_hf,4,ha);
                SetCursorLocation(g_hf,g_Scursore, g_ste);
                
                DrawCutwave(axes_cutpwave,Pdata);
                DrawCutwave(axes_cutswave,Sdata);
                
                yf = fft(Pdata(:,3),length(Pdata(:,3)));   % dataZ
                mag1 = abs(yf);
                fs = (1000/g_nsfre);
                [~,ind1] = max(mag1);
                fpw = round(ind1*fs/length(Pdata(:,3)));
                
                yff2 = fft(Sdata(:,3),length(Sdata(:,3))); % dataZ
                mag2 = abs(yff2);
                [~,ind2] = max(mag2);
                fsw = round(ind2*fs/length(Sdata(:,3)));
                
                g_PrinFreP = fpw;
                g_PrinFreS = fsw;
                
                set(hedit_p,'string',num2str(g_PrinFreP));
                set(hedit_s,'string',num2str(g_PrinFreS));
                
                DrawAllGraph();
                
            end
        end
    end

%==========================================================================
% Polarization Method Choose
%==========================================================================
% --- Executes when selected object is changed in bt_polarizmehtod.
    function btp_SelFcn(source,eventdata)
          
        switch get(eventdata.NewValue,'tag')
            case 'EVD-CM1' % Covariance matrix1
                g_polarimod = 1;
                DrawAllGraph();  
                
            case 'EVD-CM2'% Covariance matrix2
                g_polarimod = 2;
                DrawAllGraph();  
                
            case 'PCA-SVD' % SVD
                g_polarimod = 3;
                DrawAllGraph();  
                
            case 'EVD-ASM' % analytic signal
                g_polarimod = 4;
                DrawAllGraph();  
                
            otherwise
                return;
        end
           
    end
%==========================================================================
% Draw Polarization Graphs
%==========================================================================
    function DrawAllGraph()
        
        
        colorR = [0.99 0.26 0.40];
        colorB = [0 0.35 0.67];
        colorG = [0.40 0.58 0.29];
        
        
        g_ptb = GetCursorLocation(g_hf,g_Pcursorb);
        g_pte = GetCursorLocation(g_hf,g_Pcursore);
        
        g_stb = GetCursorLocation(g_hf,g_Scursorb);
        g_ste = GetCursorLocation(g_hf,g_Scursore);
        
        Pdata = g_filterdata(g_ptb:g_pte,:);
        Sdata = g_filterdata(g_stb:g_ste,:);
        
        
        cla(axes_cutpwave);
        cla(axes_cutswave);
        
        cla(axes_P3dtraj);
        cla(axes_Pazirose);
        cla(axes_PDlp);
        cla(axes_Pazimuth);
        cla(axes_Pincirose);
        cla(axes_PDpp);
        cla(axes_Pincidence);
        
        cla(axes_S3dtraj);
        cla(axes_Sazirose);
        cla(axes_SDlp);
        cla(axes_Sazimuth);
        cla(axes_Sincirose);
        cla(axes_SDpp);
        cla(axes_Sincidence);
        
        
        DrawCutwave(axes_cutpwave,Pdata);
        DrawCutwave(axes_cutswave,Sdata);
        
        delt = g_nsfre/1000; %ms-s
        ttotP = (g_pte-g_ptb+1)*delt;
        ttotS = (g_ste-g_stb+1)*delt;
        %         ttotS = (3750-3310+1)*delt;
        %--------------------------------------------------------------------------
        
        yf = fft(Pdata(:,3),length(Pdata(:,3))); % dataZ
        mag1 = abs(yf);
        fs = round(1000/g_nsfre);
        [~,ind1] = max(mag1);
%         fpw = round(ind1*fs/length(Pdata(:,3)));
        fpw = ind1*fs/length(Pdata(:,3));
        
        yff2 = fft(Sdata(:,3),length(Sdata(:,3))); % dataZ
        mag2 = abs(yff2);
        [~,ind2] = max(mag2);
%         fsw = round(ind2*fs/length(Sdata(:,3)));
        fsw = ind2*fs/length(Sdata(:,3));
        
        g_PrinFreP = fpw; % Principal frequency for P-Wave
        g_PrinFreS = fsw;
        
        set(hedit_p,'string',num2str(g_PrinFreP));
        set(hedit_s,'string',num2str(g_PrinFreS));
        
        % f = 1;  dominant frequency in Hz for earthquake
        % f = 20; % dominant frequency of synthetic time series in Hz sgy-yang
        % f = 100; % dominant frequency of synthetic time series in Hz(for myself)
        
        %--------------------------------------------------------------------------
        % Determine window to calculate polarization over, in samples
        cycs = 2; % number of cycles (2 to 3 is usually sufficient)
        if(g_PrinFreP < 0.5) % g_PrinFreP is the dominant frequency in Hz May be natural earthquake
            wndop = 200; % TW for P-wave in samples
            wndos = 200; % TW for S-wave in samples
        end
        if(0.5 <= g_PrinFreP)  % 1Hz,May be natural earthquake or low-frequency microseismic
            wndop = floor( (1/g_PrinFreP) * (1/delt) ) * cycs; % samples per cycle times # of cycles
            wndos = floor( (1/g_PrinFreS) * (1/delt) ) * cycs; % samples per cycle times # of cycles
        end
        
        %--------------------------------------------------------------------------
        twinp = (wndop)*delt;
        twins = (wndos)*delt;
        
        switch g_polarimod
            case 1 % Covariance matrix1
                [Pazi,Pinci,Pmaxeig,PDlp,PDpp]= polarize_estimation(Pdata(:,1),...
                    Pdata(:,2),Pdata(:,3),delt,ttotP,twinp);
                [Sazi,Sinci,Smaxeig,SDlp,SDpp] = polarize_estimation(Sdata(:,1),...
                    Sdata(:,2),Sdata(:,3),delt,ttotS,twins);
                
                
            case 2 % Covariance matrix2
                [Pazi,Pinci,Pmaxeig,PDlp,PDpp] = polarize_estimation2(Pdata(:,1),...
                    Pdata(:,2),Pdata(:,3),delt,ttotP,twinp);
                [Sazi,Sinci,Smaxeig,SDlp,SDpp] = polarize_estimation2(Sdata(:,1),...
                    Sdata(:,2),Sdata(:,3),delt,ttotS,twins);
                
                
            case 3 % SVD
                [Pazi,Pinci,Pmaxeig,PDlp,PDpp] = polarization_PCA(Pdata(:,1),...
                    Pdata(:,2),Pdata(:,3),delt,twinp);
                [Sazi,Sinci,Smaxeig,SDlp,SDpp] = polarization_PCA(Sdata(:,1),...
                    Sdata(:,2),Sdata(:,3),delt,twins);
                
            case 4 % analytic signal
                [Pazi,Pinci,Pmaxeig,PDlp,PDpp] = polar_analyticSig(Pdata(:,1),...
                    Pdata(:,2),Pdata(:,3),wndop);
                [Sazi,Sinci,Smaxeig,SDlp,SDpp] = polar_analyticSig(Sdata(:,1),...
                    Sdata(:,2),Sdata(:,3), wndos);
                
            otherwise
                return;
        end
        
        Draw3DPolarizTraj(axes_P3dtraj,Pdata(:,1),Pdata(:,2),Pdata(:,3),...
            'Hodogram for P-Wave');
        Draw3DPolarizTraj(axes_S3dtraj,Sdata(:,1),Sdata(:,2),Sdata(:,3),...
            'Hodogram for S-Wave');
        
  
        DrawRose(axes_Pazirose,Pazi,'Azimuth Rose for P-wave');
        DrawDlp(axes_PDlp,PDlp,'P-wave');
        DrawAzimuthHist(axes_Pazimuth,Pazi,'Azimuth histogram for P-wave');
        
        DrawRoseInc(axes_Pincirose,Pinci,'Incidence Rose for P-wave');
        DrawDpp(axes_PDpp,PDpp,'P-wave');
        DrawIncidenceHist(axes_Pincidence,Pinci,'Incidence histogram for P-wave');

        DrawRose(axes_Sazirose,Sazi,'Azimuth Rose for S-wave');
        DrawDlp(axes_SDlp,SDlp,'S-wave');
        DrawAzimuthHist(axes_Sazimuth,Sazi,'Azimuth histogram for S-wave');
        
        DrawRoseInc(axes_Sincirose,Sinci,'Incidence Rose for S-wave');
        DrawDpp(axes_SDpp,SDpp,'S-wave');
        DrawIncidenceHist(axes_Sincidence,Sinci,'Incidence histogram for S-wave');


%                 %------------------------------------------------------------------------------
%                 % Fig2 using Fig2Dlg in the paper
%                 %------------------------------------------------------------------------------
%                 colorR = [0.99 0.26 0.40];
%                 colorB = [0 0.35 0.67];
%                 colorG = [0.40 0.58 0.29];
%                 htesthandles = guihandles(Fig2Dlg);
%         
%                 cla(htesthandles.axeso);
%                 cla(htesthandles.axescp);
%                 cla(htesthandles.axescs);
%                 cla(htesthandles.axespar);
%                 cla(htesthandles.axespac);
%                 cla(htesthandles.axespas);
%                 cla(htesthandles.axessar);
%                 cla(htesthandles.axessac);
%                 cla(htesthandles.axessas);
%         
%         
%                 htest = gcf;
%                 cla(htest);
%         %         testdata = g_filterdata (11150:g_nslen,:); % for mseed data
%         %         testlen = g_nslen - 11150 + 1;
%         
%                 testdata = g_filterdata (:,:); % for synthetic data
%                 testlen = g_nslen;
%         
%                 ha = htesthandles.axeso;
%                 axes(ha);
%         
%                 plot(testdata(:,1),'color',colorR);
%                 hold on;
%                 plot(testdata(:,2),'color',colorB);
%                 hold on;
%                 plot(testdata(:,3),'color',colorG);
%                 hold off;
%                 grid on;
%                 set(ha,'ticklength',[0.005 0.0025]);
%         
%                 xlim = get(ha,'XLim');
%                 ylim = get(ha,'YLim');
%                 hxlab1 = xlabel('Samples','FontSize',14);
%                 set(hxlab1,'Position',[xlim(2)-xlim(2)*0.05,ylim(1)+ylim(2)*0.3]);
%                 hylab1 = ylabel('Amplitude','FontSize',14);
%                 set(hylab1,'Position',[xlim(1)-xlim(2)*0.015,(ylim(1)+ylim(2))*0.5]);
%         
%         
%         %         %----------------------------------
%         %         % for mseed data
%         %         %----------------------------------
%         %         Pcursorb = CreateCursor(htest,1,ha);
%         %         SetCursorLocation(htest,Pcursorb, 1000);
%         %         Pcursore = CreateCursor(htest,2,ha);
%         %         SetCursorLocation(htest,Pcursore, 1600);
%         %         Scursorb = CreateCursor(htest,3,ha);
%         %         SetCursorLocation(htest,Scursorb, 3310);
%         %         Scursore = CreateCursor(htest,4,ha);
%         %         SetCursorLocation(htest,Scursore, 3750);
%         %         Pdata = testdata(1000:1600,:);
%         %         Sdata = testdata(3310:3750,:);
%         
%                 %----------------------------------
%                 % for sythetic data
%                 %----------------------------------
%                 Pcursorb = CreateCursor(htest,1,ha);
%                 SetCursorLocation(htest,Pcursorb, g_ptb);
%                 Pcursore = CreateCursor(htest,2,ha);
%                 SetCursorLocation(htest,Pcursore, g_pte);
%                 Scursorb = CreateCursor(htest,3,ha);
%                 SetCursorLocation(htest,Scursorb, g_stb);
%                 Scursore = CreateCursor(htest,4,ha);
%                 SetCursorLocation(htest,Scursore, g_ste);
%         
%                 Pdata = testdata(g_ptb:g_pte,:);
%                 Sdata = testdata(g_stb:g_ste,:);
%                 
%                 
%                 colorR2 = [0.58 0.16 0.13];
% 
%                 ha = htesthandles.axescp;
%                 axes(ha);
%         
%                 plot(Pdata(:,1),'color',colorR,'LineWidth',2);
%                 hold on;
%                 plot(Pdata(:,2),'color',colorB,'LineWidth',2);
%                 hold on;
%                 plot(Pdata(:,3),'color',colorG,'LineWidth',2);
%                 hold off;
%                 grid on;
%                 set(ha,'ticklength',[0.005 0.0025]);
%         
%                 xlim = get(ha,'XLim');
%                 ylim = get(ha,'YLim');
%                 hxlab1 = xlabel('Samples','FontSize',14);
%                 set(hxlab1,'Position',[xlim(2)-xlim(2)*0.05,ylim(1)+ylim(2)*0.3]);
%                 hylab1 = ylabel('Amplitude','FontSize',14);
%                 set(hylab1,'Position',[xlim(1)-xlim(2)*0.015,(ylim(1)+ylim(2))*0.5]);
%                 refreshdata;
%                 drawnow;
%                 
%                 
%                 ha = htesthandles.axescs;
%                 axes(ha);
%                 
%                 plot(Sdata(:,1),'color',colorR,'LineWidth',2);
%                 hold on;
%                 plot(Sdata(:,2),'color',colorB,'LineWidth',2);
%                 hold on;
%                 plot(Sdata(:,3),'color',colorG,'LineWidth',2);
%                 hold off;
%                 grid on;
%                 set(ha,'ticklength',[0.005 0.0025]);
%         
%                 xlim = get(ha,'XLim');
%                 ylim = get(ha,'YLim');
%                 hxlab1 = xlabel('Samples','FontSize',14);
%                 set(hxlab1,'Position',[xlim(2)-xlim(2)*0.05,ylim(1)+ylim(2)*0.3]);
%                 hylab1 = ylabel('Amplitude','FontSize',14);
%                 set(hylab1,'Position',[xlim(1)-xlim(2)*0.015,(ylim(1)+ylim(2))*0.5]);
%                 refreshdata;
%                 drawnow;
%         
%                 
%                 ha = htesthandles.axespar;
%                 axes(ha);
%                 polar_plot_own2(Pazi);
%                 title('Azimuth rose for P-wave','color',colorB);
%         
%                 ha = htesthandles.axespac;
%                 axes(ha);
%                 histogram(Pazi,'LineWidth',0.5,'FaceColor',colorB,'EdgeColor',colorR2); hold on;
%                 xlabel('Azimuth / Degree','FontSize',12);ylabel('Numbers','FontSize',12);
%                 title('Azimuth histogram for P-wave','color',colorB);
%         
%                 ha = htesthandles.axespas;
%                 axes(ha);
%                 h1 = plot(Pazi,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%                 set(h1,'Marker','o');
%                 set(h1,'MarkerEdgeColor',colorR);
%                 set(h1,'MarkerFaceColor',colorB);
%                 xlabel('Nums of TW');ylabel('Degree');
%                 title('Azimuth statistic for P-wave','color',colorB);
%                 grid on;
%         
%                 ha = htesthandles.axessar;
%                 axes(ha);
%                 polar_plot_own2(Sazi);
%                 title('Azimuth rose for S-wave','color',colorB);
%         
%                 ha = htesthandles.axessac;
%                 axes(ha);
%                 histogram(Sazi,'LineWidth',0.5,'FaceColor',colorB,'EdgeColor',colorR2); hold on;
%                 xlabel('Azimuth / Degree','FontSize',12);ylabel('Numbers','FontSize',12);
%                 title('Azimuth histogram for S-wave','color',colorB);
%         
%                 ha = htesthandles.axessas;
%                 axes(ha);
%                 h2 = plot(Sazi,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%                 set(h2,'Marker','o');
%                 set(h2,'MarkerEdgeColor',colorR);
%                 set(h2,'MarkerFaceColor',colorB);
%                 xlabel('Nums of TW');ylabel('Degree');
%                 title('Azimuth statistic for S-wave','color',colorB);
%                 grid on;
%         
%                 refreshdata;
%                 drawnow;
%         
%                 %------------------------------------------------------------------------------
%                 % Fig2- END
%                 %------------------------------------------------------------------------------
%         
%                 %------------------------------------------------------------------------------
%                 % Fig3-Fig4
%                 % P-Wave and S-wave Azimuth comparison
%                 %------------------------------------------------------------------------------
%         
%                 figure(2);
%                 set(gcf,'color','w');
%                 subplot(2, 3, 1);
%                 polar_plot_own2(Pazi);
%                 title('Azimuth rose for P-wave','color',colorB);
%         
%                 subplot(2, 3, 2);
%                 histogram(Pazi,'LineWidth',0.5,'FaceColor',colorB,'EdgeColor',colorR2); hold on;
%                 xlabel('Azimuth / Degree','FontSize',12);ylabel('Numbers','FontSize',12);
%                 title('Azimuth histogram for P-wave','color',colorB);
%         
%                 subplot(2, 3, 3);
%                 h1 = plot(Pazi,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%                 set(h1,'Marker','o');
%                 set(h1,'MarkerEdgeColor',colorR);
%                 set(h1,'MarkerFaceColor',colorB);
%                 xlabel('Nums of TW');ylabel('Degree');
%                 title('Azimuth statistic for P-wave','color',colorB);
%                 grid on;
%         
%                 subplot(2, 3, 4);
%                 polar_plot_own2(Sazi);
%                 title('Azimuth rose for S-wave','color',colorB);
%         
%                 subplot(2, 3, 5);
%                 histogram(Sazi,'LineWidth',0.5,'FaceColor',colorB,'EdgeColor',colorR2); hold on;
%                 xlabel('Azimuth / Degree','FontSize',12);ylabel('Numbers','FontSize',12);
%                 title('Azimuth histogram for S-wave','color',colorB);
%         
%                 subplot(2, 3, 6);
%                 h2 = plot(Sazi,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%                 set(h2,'Marker','o');
%                 set(h2,'MarkerEdgeColor',colorR);
%                 set(h2,'MarkerFaceColor',colorB);
%                 xlabel('Nums of TW');ylabel('Degree');
%                 title('Azimuth statistic for S-wave','color',colorB);
%                 grid on;
%                 %------------------------------------------------------------------------------
%                 % Fig3-Fig4 END
%                 %------------------------------------------------------------------------------
%         
%         
%         %-----------------------
%         % Time consumption
%         tic;
%         [Pazi1,Pinci1,Pmaxeig1,PDlp1,PDpp1] = polarize_estimation(Pdata(:,1),...
%             Pdata(:,2),Pdata(:,3),delt,ttotP,twinp);
%         tCCM1 = toc;
%         
%         tic;
%         [Pazi2,Pinci2,Pmaxeig2,PDlp2,PDpp2] = polarize_estimation2(Pdata(:,1),...
%             Pdata(:,2),Pdata(:,3),delt,ttotP,twinp);
%         tCCM2 = toc;
%         
%         tic;
%         [Pazi3,Pinci3,Pmaxeig3,PDlp3,PDpp3] = polarization_PCA(Pdata(:,1),...
%             Pdata(:,2),Pdata(:,3),delt,twinp);
%         tSVD = toc;
%         
%         tic;
%         [Pazi4,Pinci4,Pmaxeig4,PDlp4,PDpp4] = polar_analyticSig(Pdata(:,1),...
%             Pdata(:,2),Pdata(:,3),wndop);
%         tCM = toc;
%         
%         fprintf('\tTime consumption:\n');
%         fprintf('\ttCCM1\t tCCM2\t tSVD\t tAnaSig\n');
%         fprintf('\t%f\t %f\t %f\t %f\n',tCCM1,tCCM2,tSVD,tCM);
%         fprintf('\t！！！！！！！！！！！！！！！！！！！！！！！！\n');
%         
%         max_Pazi1 = max(Pazi1); %max
%         min_Pazi1 = min(Pazi1); %min
%         diff1 = max_Pazi1 - min_Pazi1; %max-min
%         ave_Pazi1 = mean(Pazi1); %average
%         
%         max_Pazi2 = max(Pazi2); %max
%         min_Pazi2 = min(Pazi2); %min
%         diff2 = max_Pazi2 - min_Pazi2; %max-min
%         ave_Pazi2 = mean(Pazi2); %average
%         
%         max_Pazi3 = max(Pazi3); %max
%         min_Pazi3 = min(Pazi3); %min
%         diff3 = max_Pazi3 - min_Pazi3; %max-min
%         ave_Pazi3 = mean(Pazi3); %average
%         
%         max_Pazi4 = max(Pazi4); %max
%         min_Pazi4 = min(Pazi4); %min
%         diff4 = max_Pazi4 - min_Pazi4; %max-min
%         ave_Pazi4 = mean(Pazi4); %average
%         
%         fprintf('\tAzi Max value:\n');
%         fprintf('\tCCM1\t\t CCM2\t\t SVD\t\t AnaSig\n');
%         fprintf('\t%.20f\t %.20f\t %.20f\t %.20f\n',max_Pazi1,max_Pazi2,max_Pazi3,max_Pazi4);
%         fprintf('\t！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！\n');
%         
%         fprintf('\tAzi Min value:\n');
%         fprintf('\tCCM1\t\t CCM2\t\t SVD\t\t AnaSig\n');
%         fprintf('\t%.20f\t %.20f\t %.20f\t %.20f\n',min_Pazi1,min_Pazi2,min_Pazi3,min_Pazi4);
%         fprintf('\t！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！\n');
%         
%         fprintf('\tAzi Diff value for max-min:\n');
%         fprintf('\tCCM1\t\t CCM2\t\t SVD\t\t AnaSig\n');
%         fprintf('\t%.20f\t %.20f\t %.20f\t %.20f\n',diff1,diff2,diff3,diff4);
%         fprintf('\t！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！\n');
%         
%         fprintf('\tAzi average value:\n');
%         fprintf('\tCCM1\t\t CCM2\t\t SVD\t\t AnaSig\n');
%         fprintf('\t%.20f\t %.20f\t %.20f\t %.20f\n',ave_Pazi1,ave_Pazi2,ave_Pazi3,ave_Pazi4);
%         fprintf('\t！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！\n');
% 		
%         fprintf('\tlength of the calculated data:%d\n\n',size(Pdata,1));
%         
%         
%         %-----------------------
%         %-----------------------
%         % Azimuth Comparison using rose for Fig.6
%         figure(3);
%         set(gcf,'color','w');
%         subplot(2, 2, 1);
%         polar_plot_own2(Pazi1);
%         title('Azimuth rose for EVD-CM1','color',colorB);
%         
%         subplot(2, 2, 2);
%         polar_plot_own2(Pazi2);
%         title('Azimuth rose for EVD-CM2','color',colorB);
%         
%         subplot(2, 2, 3);
%         polar_plot_own2(Pazi3);
%         title('Azimuth rose for PCA-SVD','color',colorB);
%         
%         subplot(2, 2, 4);
%         polar_plot_own2(Pazi4);
%         title('Azimuth rose for EVD-ASM','color',colorB);
%         
%         
%         %-----------------------
%         % Incidence Comparison using rose for Fig.7
%         
%         figure(4);
%         set(gcf,'color','w');
%         subplot(2, 2, 1);
%         polar_plot_In_own2(Pinci1);
%         title('Incidence rose for EVD-CM1','color',colorB);
%         
%         subplot(2, 2, 2);
%         polar_plot_In_own2(Pinci2);
%         title('Incidence rose for EVD-CM2','color',colorB);
%         
%         subplot(2, 2, 3);
%         polar_plot_In_own2(Pinci3);
%         title('Incidence rose for PCA-SVD','color',colorB);
%         
%         subplot(2, 2, 4);
%         polar_plot_In_own2(Pinci4);
%         title('Incidence rose for EVD-ASM','color',colorB);
%         
%         %-----------------------
%         % Azimuth Comparison using statistical graph for Fig.7
%         
%         figure(5);
%         yMin = min(min([Pazi1;Pazi2;Pazi3;Pazi4]));yMax = max(max([Pazi1;Pazi2;Pazi3;Pazi4]));
%         subplot(2, 2, 1);
%         h2 = plot(Pazi1,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[yMin yMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree');
%         title('Azimuth statistic for EVD-CM1','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 2);
%         h2 = plot(Pazi2,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[yMin yMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree');
%         title('Azimuth statistic for EVD-CM2','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 3);
%         h2 = plot(Pazi3,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[yMin yMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree');
%         title('Azimuth statistic for PCA-SVD','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 4);
%         h2 = plot(Pazi4,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[yMin yMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree');
%         title('Azimuth statistic for EVD-ASM','color',colorB);
%         grid on;
%         
%         
%         %-----------------------
%         % Azimuth histogram graph for Fig.7
%         figure(6);
%         colorR2 = [0.58 0.16 0.13];
%         subplot(2, 2, 1);
%         histogram(Pazi1,'LineWidth',0.5,'FaceColor',colorB,'EdgeColor',colorR2); hold on;
%         xlabel('Azimuth / Degree','FontSize',12);ylabel('Numbers','FontSize',12);
%         title('Azimuth histogram for EVD-CM1','color',colorB);
%  
%         grid on;
% 
%         subplot(2, 2, 2);
%         histogram(Pazi2,'LineWidth',0.5,'FaceColor',colorB,'EdgeColor',colorR2); hold on;
%         xlabel('Azimuth / Degree','FontSize',12);ylabel('Numbers','FontSize',12);
%         title('Azimuth histogram for EVD-CM2','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 3);
%         histogram(Pazi3,'LineWidth',0.5,'FaceColor',colorB,'EdgeColor',colorR2); hold on;
%         xlabel('Azimuth / Degree','FontSize',12);ylabel('Numbers','FontSize',12);
%         title('Azimuth histogram for PCA-SVD','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 4);
%         histogram(Pazi4,'LineWidth',0.5,'FaceColor',colorB,'EdgeColor',colorR2); hold on;
%         xlabel('Azimuth / Degree','FontSize',12);ylabel('Numbers','FontSize',12);
%         title('Azimuth histogram for EVD-ASM','color',colorB);
%         grid on;
%         
%         
%         %-----------------------
%         % Degree of rectilinearity Comparison using statistical graph
%         figure(7);
%         dlpMin = min(min([PDlp1;PDlp2;PDlp3;PDlp4]));dlpMax = max(max([PDlp1;PDlp2;PDlp3;PDlp4]));
%         subplot(2, 2, 1);
%         h2 = plot(PDlp1,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[dlpMin dlpMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree of rectilinearity');
%         title('Dlp statistic for EVD-CM1','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 2);
%         h2 = plot(PDlp2,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[dlpMin dlpMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree of rectilinearity');
%         title('Dlp statistic for EVD-CM2','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 3);
%         h2 = plot(PDlp3,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[dlpMin dlpMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree of rectilinearity');
%         title('Dlp statistic for PCA-SVD','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 4);
%         h2 = plot(PDlp4,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[dlpMin dlpMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree of rectilinearity');
%         title('Dlp statistic for EVD-ASM','color',colorB);
%         grid on;
%         
%         
%          %-----------------------
%         %  The degree of planarity Comparison using statistical graph
%         figure(8);
%         dppMin = min(min([PDpp1;PDpp2;PDpp3;PDpp4]));dppMax = max(max([PDpp1;PDpp2;PDpp3;PDpp4]));
%         subplot(2, 2, 1);
%         h2 = plot(PDpp1,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[dppMin dppMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree of planarity');
%         title('Dpp statistic for EVD-CM1','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 2);
%         h2 = plot(PDpp2,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[dppMin dppMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree of planarity');
%         title('Dpp statistic for EVD-CM2','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 3);
%         h2 = plot(PDpp3,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[dppMin dppMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree of planarity');
%         title('Dpp statistic for PCA-SVD','color',colorB);
%         grid on;
%         
%         subplot(2, 2, 4);
%         h2 = plot(PDpp4,'-','LineWidth',1,'color',colorB,'MarkerSize',2); hold on;
%         set(gca,'YLim',[dppMin dppMax]);
%         set(h2,'Marker','o');
%         set(h2,'MarkerEdgeColor',colorR);
%         set(h2,'MarkerFaceColor',colorB);
%         xlabel('Nums of TW');ylabel('Degree of planarity');
%         title('Dpp statistic for EVD-ASM','color',colorB);
%         grid on;
% 
%         %------------------------------------------------------------------------------
%         % Fig1-hodogram using testDlg in the paper
%         %------------------------------------------------------------------------------
%                 colorR = [0.99 0.26 0.40];
%                 colorB = [0 0.35 0.67];
%                 colorG = [0.40 0.58 0.29];
%                 htestdlghandles = guihandles(TestDlg);
%         
%                 ha = htestdlghandles.axes_owave;
%         
%                 htest = gcf;
%         %         testdata = g_filterdata (11150:g_nslen,:); % for real data(.sac)
%         %         testlen = g_nslen - 11150 + 1;
%         
%                 testdata = g_filterdata;
%                 axes(ha);
%         
%                 plot(testdata(:,1),'color',colorR);
%                 hold on;
%                 plot(testdata(:,2),'color',colorB);
%                 hold on;
%                 plot(testdata(:,3),'color',colorG);
%                 hold off;
%                 grid on;
%                 set(ha,'ticklength',[0.005 0.0025]);
%         
%                 xlim = get(ha,'XLim');
%                 ylim = get(ha,'YLim');
%                 hxlab1 = xlabel('Samples','FontSize',10);
%                 set(hxlab1,'Position',[xlim(2)-xlim(2)*0.05,ylim(1)+ylim(2)*0.3]);
%                 hylab1 = ylabel('Amplitude','FontSize',10);
%                 set(hylab1,'Position',[xlim(1)-xlim(2)*0.015,(ylim(1)+ylim(2))*0.5]);
%         
%                 refreshdata;
%                 drawnow;
%         
%                 %----------------------------------
%         
%         
%         
%                 Pcursorb = CreateCursor(htest,1,ha);
%                 SetCursorLocation(htest,Pcursorb, g_ptb);
%                 Pcursore = CreateCursor(htest,2,ha);
%                 SetCursorLocation(htest,Pcursore, g_pte);
%         
%                 Pdata = testdata(g_ptb:g_pte,:);
%         
%                 ha = htestdlghandles.axes_cutwave;
%                 axes(ha);
%         
%                 plot(Pdata(:,1),'color',colorR,'LineWidth',2);
%                 hold on;
%                 plot(Pdata(:,2),'color',colorB,'LineWidth',2);
%                 hold on;
%                 plot(Pdata(:,3),'color',colorG,'LineWidth',2);
%                 hold off;
%                 grid on;
%                 set(ha,'ticklength',[0.005 0.0025]);
%         
%                 xlim = get(ha,'XLim');
%                 ylim = get(ha,'YLim');
%                 hxlab1 = xlabel('Samples','FontSize',10);
%                 set(hxlab1,'Position',[xlim(2)-xlim(2)*0.05,ylim(1)+ylim(2)*0.3]);
%                 hylab1 = ylabel('Amplitude','FontSize',10);
%                 set(hylab1,'Position',[xlim(1)-xlim(2)*0.015,(ylim(1)+ylim(2))*0.5]);
%                 refreshdata;
%                 drawnow;
%         
%                 %----------------------------------
%                 % 3D trajectory
%                 %----------------------------------
%         
%                 ptrace = Pdata;
%                 pfact = max(sqrt(sum(ptrace'.^2,2))); %Normalization
%                 ptrace_norm = ptrace/pfact;
%                 ha = htestdlghandles.axes_3d;
%                 axes(ha);
%                 plot_dir3N (ptrace_norm(:,1),ptrace_norm(:,2),ptrace_norm(:,3));
%                 % plot_dir3 (vp(1,:)',vp(2,:)',vp(3,:)');
%         
%                 meanX = mean(ptrace_norm,1);
%                 [coeff,score,latent] = pca(ptrace_norm,'Algorithm','SVD');
%                 EigVect = coeff(:,1);
%                 tAZ = [min(score(:,1))-.2, max(score(:,1))+.2]; 
%                 endptsp = [meanX + tAZ(1)*EigVect'; meanX + tAZ(2)*EigVect']; % Recovered PCA data
%                 hold on;
%         
%                 plot3 (endptsp(:,1),endptsp(:,2),endptsp(:,3),'--','Color',colorR,'LineWidth',3);
%                 xlabel('X');ylabel('Y');zlabel('Z');
%         
%                 grid on;
%                 hold off;
%         
%                 %----------------------------------
%                 % 2D trajectory
%                 %----------------------------------
%         
%                 ptrace = Pdata;
%                 pfact = max(sqrt(sum(ptrace'.^2,2))); %Normalization
%                 ptrace_norm = ptrace/pfact;
%                 ha = htestdlghandles.axes_2d;
%                 axes(ha);
%                 h1 = plot(ptrace_norm(:,1),ptrace_norm(:,2),'-','LineWidth',1,'Color',colorB,'MarkerSize',3);
%                 hold on;
%                 set(h1,'Marker','.');
%                 set(h1,'MarkerEdgeColor','r');
%                 set(h1,'MarkerFaceColor',colorB);
%         
%                 meanX = mean(ptrace_norm,1);
%                 [coeff,score,latent] = pca(ptrace_norm,'Algorithm','SVD');
%                 EigVect = coeff(:,1);
%                 tAZ = [min(score(:,1))-.2, max(score(:,1))+.2];
%                 endptsp = [meanX + tAZ(1)*EigVect'; meanX + tAZ(2)*EigVect']; % Recovered PCA data
%                 hold on;
%         
%                 plot (endptsp(:,1),endptsp(:,2),'--','Color',colorR,'LineWidth',3);
%                 xlabel('X');ylabel('Y');
%         
%                 grid on;
%                 hold off;
%         %------------------------------------------------------------------------------
%         % Hodogram END
%         %------------------------------------------------------------------------------
     
    end

%==========================================================================
%==========================================================================
    function Drawwave(ha,alldata)
        
        
        colorR = [0.86 0.30 0.43];
        % colorB = [0.18 0.66 0.87];
        % colorG = [0.40 0.58 0.29];
        
        axes(ha);
        
        plot(alldata(:,1),'color',colorR,'LineWidth',1);
        hold on;
        plot(alldata(:,2),'w','LineWidth',1);
        hold on;
        hpltz = plot(alldata(:,3),'c','LineWidth',1);
        hold off;
        set(ha,'color', g_backColor, 'Ycolor','w', 'Xcolor','w');
        
        set(ha,'ticklength',[0.005 0.0025]);
        
        xlim = get(ha,'XLim');
        ylim = get(ha,'YLim');
        hxlab1 = xlabel('Samples','FontSize',12,'Color','w');
        set(hxlab1,'Position',[xlim(2)-xlim(2)*0.05,ylim(1)+ylim(2)*0.3]);
        
        hylab1 = ylabel('Amplitude','FontSize',12,'Color','w');
        set(hylab1,'Position',[xlim(1)-xlim(2)*0.02,(ylim(1)+ylim(2))*0.5]);
        grid on;
        refreshdata;
        drawnow;
    end

%==========================================================================
%==========================================================================
    function DrawCutwave(ha,alldata)
        
        
        colorR = [0.86 0.30 0.43];
        % colorB = [0.18 0.66 0.87];
        % colorG = [0.40 0.58 0.29];
        
        axes(ha);
        
        plot(alldata(:,1),'color',colorR,'LineWidth',1);
        hold on;
        plot(alldata(:,2),'w','LineWidth',1);
        hold on;
        plot(alldata(:,3),'c','LineWidth',1);
        hold off;
        set(ha,'color', g_backColor, 'Ycolor','w', 'Xcolor','w');
        
        set(ha,'ticklength',[0.005 0.0025]);
        
        xlim = get(ha,'XLim');
        ylim = get(ha,'YLim');
        hxlab1 = xlabel('Samples','FontSize',12,'Color','w');
        set(hxlab1,'Position',[xlim(2)-xlim(2)*0.05,ylim(1)+ylim(2)*0.3]);
        
        hylab1 = ylabel('Amplitude','FontSize',12,'Color','w');
        set(hylab1,'Position',[xlim(1)-xlim(2)*0.05,(ylim(1)+ylim(2))*0.5]);
        grid on;
        refreshdata;
        drawnow;
    end
%==========================================================================
    function Draw3DPolarizTraj(ha,X,Y,Z,titlename)
        ptrace = [X Y Z];
        pfact = max(sqrt(sum(ptrace'.^2,2))); %Normalization
        ptrace_norm = ptrace/pfact;
        axes(ha);
        plot_dir3 (ptrace_norm(:,1),ptrace_norm(:,2),ptrace_norm(:,3));
        % plot_dir3 (vp(1,:)',vp(2,:)',vp(3,:)');
        
        meanX = mean(ptrace_norm,1);
        [coeff,score,latent] = pca(ptrace_norm,'Algorithm','SVD');
        EigVect = coeff(:,1);
        tAZ = [min(score(:,1))-.2, max(score(:,1))+.2];
        endptsp = [meanX + tAZ(1)*EigVect'; meanX + tAZ(2)*EigVect']; % Recovered PCA data
        hold on;
        plot3 (endptsp(:,1),endptsp(:,2),endptsp(:,3),'y--','LineWidth',3);
        xlabel('E');ylabel('N');zlabel('Z');
        title(titlename,'color','w');
        grid on;
        hold off;
        set(ha,'color', g_backColor, 'Ycolor','w', 'Xcolor','w', 'Zcolor','w');

    end
%=======================================================================
% Azimuth/Incidence histogram
%=======================================================================
    function DrawAzimuthHist(ha,azi,titlename)
        
        axes(ha);
        histogram(azi,'LineWidth',0.5,'FaceColor','c','EdgeColor','y'); hold on;
        xlabel('Azimuth / Degree','color','w','FontSize',12);ylabel('Numbers','color','w','FontSize',12);
        title(titlename,'color','w');
        set(ha,'color', g_backColor, 'Ycolor','w', 'Xcolor','w');
        grid on;
    end

    function DrawIncidenceHist(ha,inci,titlename)
        axes(ha);
        histogram(inci,'LineWidth',0.5,'FaceColor','c','EdgeColor','y'); hold on;
        xlabel('Incidence / Degree','color','w','FontSize',12);ylabel('Numbers','color','w','FontSize',12);
        title(titlename,'color','w');
        set(ha,'color', g_backColor, 'Ycolor','w', 'Xcolor','w');
        grid on;
    end

%=======================================================================
% degree of linear polarization�� The degree of planarity
%=======================================================================
    function DrawDlp(ha,Dlp,titlename)
        
        axes(ha);
        h1 = plot(Dlp,'c-.','LineWidth',0.5,'MarkerSize',2); hold on;
        set(h1,'Marker','o');
        set(h1,'MarkerEdgeColor','y');
        set(h1,'MarkerFaceColor','r');
        xlabel('Nums of TW','color','w','FontSize',12);ylabel('Degree of rectilinearity','color','w','FontSize',12);
        title(titlename,'color','w');
        set(ha,'color', g_backColor, 'Ycolor','w', 'Xcolor','w');
        grid on;
    end

    function DrawDpp(ha,Dpp,titlename)
        axes(ha);
        h1 = plot(Dpp,'c-.','LineWidth',0.5,'MarkerSize',2); hold on;
        set(h1,'Marker','o');
        set(h1,'MarkerEdgeColor','y');
        set(h1,'MarkerFaceColor','r');
        xlabel('Nums of TW','color','w','FontSize',12);ylabel('Degree of planarity','color','w','FontSize',12);
        title(titlename,'color','w');
        grid on;
        set(ha,'color', g_backColor, 'Ycolor','w', 'Xcolor','w');
    end
%=======================================================================
    function DrawRose(ha,azi,titlename)
        polar_plot_own(ha,azi);
        
        set(ha,'color',g_backColor);
        set(ha,'xcolor',g_backColor,'ycolor',g_backColor);
        ha.Title.String = titlename;
        ha.Title.Color = 'w';
    end
%=======================================================================
    function DrawRoseInc(ha,azi,titlename)
        polar_plot_In_own(ha,azi);
        set(ha,'color',g_backColor);
        set(ha,'xcolor',g_backColor,'ycolor',g_backColor);
        ha.Title.String = titlename;
        ha.Title.Color = 'w';
    end
%=======================================================================
%     function DrawScatter(ha,dp,titlename)
%         
%         udp = unique(dp*pi/180);
%         for i=1:length(udp)
%             dpnum(i)=length(find((dp*pi/180)==udp(i))); %拷窃由柴
%         end
%         psize = 20;
%         polarscatter(ha,udp,dpnum,psize,'filled','c');
%         pax = ha;
%         set(pax,'color',g_backColor);
%         pax.LineWidth = 1;
%         pax.GridColor = 'w';
%         pax.ThetaColor = 'w';
%         pax.RColor = 'w';
% %         title(titlename,'color','w');
%         pax.Title.String = titlename;
%         pax.Title.Color = 'w';
%     end
    
%=======================================================================
    function pb_r3d_Callback()
        hRot = rotate3d;
        if(IsOpenR3D == 0)
            hRot.Enable = 'on';
            set(hbtn_r3d,'string','R3DOn');
            IsOpenR3D = 1;
            setAllowAxesRotate(hRot,axes_swave,false); % disable rotating for other plot
            setAllowAxesRotate(hRot,axes_cutpwave,false);
            setAllowAxesRotate(hRot,axes_cutswave,false);
            setAllowAxesRotate(hRot,axes_PDlp,false);
            setAllowAxesRotate(hRot,axes_Pazimuth,false);
            setAllowAxesRotate(hRot,axes_PDpp,false);
            setAllowAxesRotate(hRot,axes_Pincidence,false);
            setAllowAxesRotate(hRot,axes_SDlp,false);
            setAllowAxesRotate(hRot,axes_Sazimuth,false);
            setAllowAxesRotate(hRot,axes_SDpp,false);
            setAllowAxesRotate(hRot,axes_Sincidence,false);
            
        elseif(IsOpenR3D == 1)
            hRot.Enable = 'off';
            set(hbtn_r3d,'string','R3DOff');
            IsOpenR3D = 0;         
        end 
    end

end
%=======================================================================
%=======================================================================
