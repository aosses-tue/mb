function VoD_stage_4_to_5(directory)
% function VoD_stage_4_to_5(directory)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file name: r20150220_update.m
% Created on    : 23/02/2015
% Last update on: 23/02/2015 % Update this date manually
% Last use on   : 23/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 1;
Diary(mfilename,bDiary);

if nargin < 1
    %directory = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-sec-Harm\'; % including Harmonic
    directory = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-sec-Harm-G1\'; % including Harmonic
end

dir_src_model = [directory 'Stage4' delim];
dir_dst_model = [directory 'Stage5' delim];

Mkdir(dir_dst_model);

% dir_src_meas  = [root_dir '02-Wav-files\2015-02-wav-files\01-as-obtained\'];
% dir_dst_meas  = [root_dir '02-Wav-files\2015-02-wav-files\02-calibrated\'  ];

% Original name         Destination name
% mode-1-v_1-ane.wav    model-ac-2-close-ane.wav    where model/meas; 1 = ac mode 2; 1 = close
% mode-3-v_2-rev.wav    meas-ac-4-dist-rev.wav      where model/meas; 3 = ac mode 4; 2 = distant

callvl = 100; % 100 = AMT, 90 = Fastl
acmode = [2 4];
lvlmode= [60 78]; 
mic    = [1 2]; % far - close
type   = {'ane'; 'rev'}; 
typefull = {'anechoic';'reverberant'};

for i = 1:length(acmode)

    % 1. Obtain level differences deltaCloseFar
    % 2. Calibrate and rename model signals
    for j = 1:length(type)
        
        callvl_dBFS = lvlmode(i)-callvl;
        
        fname1 = sprintf('%smode-%.0f-v_1-%s.wav',dir_src_model,acmode(i)-1,type{j});
        fname2 = sprintf('%smode-%.0f-v_2-%s.wav',dir_src_model,acmode(i)-1,type{j});
        
        [x1 fs] = Wavread(fname1);
        rmsFar = rmsdb(x1);
        
        [x2 fs] = Wavread(fname2);
        rmsClose = rmsdb(x2);
        deltaCloseFar(j,i) = rmsClose - rmsFar;
        
        corr(i) = callvl_dBFS - rmsClose;
        
        x1 = setdbspl(x1,rmsClose+100);
        y1 = From_dB( corr(i)-deltaCloseFar(j,i) )*x1;
        y2 = From_dB( corr(i)                    )*x2;
        
        fname1o = sprintf('%smodel-ac-%.0f-dist-%s.wav',dir_dst_model,acmode(i),type{j});
        fname2o = sprintf('%smodel-ac-%.0f-close-%s.wav' ,dir_dst_model,acmode(i),type{j});
        
        Wavwrite(y1,fs,fname1o);
        Wavwrite(y2,fs,fname2o);
        
    end
    
    % for j = 1:length(type)
    % 
    %     callvl_dBFS = lvlmode(i)-callvl;
    % 
    %     fname1 = sprintf('%smeas-ac-mode-%.0f-close-%s.wav',dir_src_meas,acmode(i),typefull{j}); % meas-ac-mode-2-close-anechoic.wav
    %     fname2 = sprintf('%smeas-ac-mode-%.0f-dist-%s-HP.wav',dir_src_meas,acmode(i),typefull{j});
    % 
    %     [x1 fs] = Wavread(fname1);
    % 
    %     [x2 fs] = Wavread(fname2);
    %     x2 = setdbspl(x2,callvl_dBFS+100);
    % 
    %     x1 = setdbspl(x1,rmsdb(x2)+100);
    %     y1 = From_dB( -deltaCloseFar(j,i) )*x1;
    %     y2 =                               x2;
    % 
    %     fname1o = sprintf('%smeas-ac-%.0f-dist-%s.wav',dir_dst_meas,acmode(i),type{j});
    %     fname2o = sprintf('%smeas-ac-%.0f-close-%s.wav' ,dir_dst_meas,acmode(i),type{j});
    % 
    %     Wavwrite(y1,fs,fname1o);
    %     Wavwrite(y2,fs,fname2o);
    % 
    % end
    
end

disp('')

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
