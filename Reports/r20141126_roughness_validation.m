function outs = r20141126_roughness_validation(options)
% function outs = r20141126_roughness_validation(options)
%
% 1. Description:
%       Implement and validate the model of Roughness (off-line).
%       Run this script first setting bCreate to 1
%       Calibration: RMS of -30 dBFS corresponds to an average level of 60 dB
% 
% 2. Stand-alone example:
%       options.bDiary = 0; % To generate log-file
%       options.bCreate = 0; % Set to 1 if you want to generate the test signals
%       options.bDoExp0 = 1; % Fastl2007, Fig. 11.1. Daniel1997, Fig.5. 
%                            % (Only the first 8192 samples are used)
%       r20141126_roughness_validation(options);
%       
%       options = [];
%       options.bCreate = 0;
%       outs = r20141126_roughness_validation(options);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 25/11/2014
% Last update on: 25/06/2015 % Update this date manually
% Last use on   : 14/07/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    options = [];
end
options = ef(options,'bDiary', 0);

options = ef(options,'bDoExp0',1);
options = ef(options,'bDoExp1',0);
options = ef(options,'bDoExp2',0);
options = ef(options,'bDoExp3',0);
options = ef(options,'bDoExp4',0);
options = ef(options,'bDoExp5',0);
options = ef(options,'bDoExp6',0);

bDiary = options.bDiary;
Diary(mfilename,bDiary);

close all

pathaudio_src = [Get_TUe_data_paths('db_Fastl2007')];
pathaudio   = [Get_TUe_paths('outputs') 'Fastl2007_R_test_20141126' delim];
pathaudio_D = [Get_TUe_paths('outputs') 'Daniel1997_R_test_20141126' delim];

p = Get_date;
pathfigures = [Get_TUe_paths('outputs') 'Figures-' p.date4files delim];
Mkdir(pathfigures);
listed_files = {};
count_files = 1;

h = [];

N           = 8192;
optsDebug.all = 1;
optsDebug.ki  = 1;
bCreate     = options.bCreate;

if bCreate
    
    opts.bDoRough = 1;
    opts.dur = (20*N/44100); %  approx. 200e-3;
    opts.bGen_test_tones = 1;
    
    opts.bPsySound = 1;
    opts.bDoRamp = 0; 
    opts.bDoZeroPadding = 0;
    
    outs = Generate_reference_sounds_Zwicker2007(opts);
    outs2 = Generate_reference_sounds(opts);
    
end

bDoExp0 = options.bDoExp0; % Reference
bDoExp1 = options.bDoExp1; % R(m)         - Fastl2007, Fig 11.1. Daniel1997, Fig 5
bDoExp2 = options.bDoExp2; % R(fc,fmod)   - Fastl2007, Fig 11.2. Daniel1997, Fig 3
bDoExp3 = options.bDoExp3; % R(type,fmod) - Fastl2007, Fig 11.3
bDoExp4 = options.bDoExp4; % R(type,SPL)  - Fastl2007, Fig 11.4
bDoExp5 = options.bDoExp5; % R(fdev)      - Fastl2007, Fig.11.5
bDoExp6 = options.bDoExp6; % R(fmod), FM  - Daniel1997, Fig 9

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenames0  = {'rough_ref'};

filenames1  = {'track_38_t01', 1;...
               'track_38_t02', 0.7; ...
               'track_38_t03', 0.4; ...
               'track_38_t04', 0.25; ...
               'track_38_t05', 0.125; ...
               'track_38_t06', 0.1; ...
               'track_38_t07', 0};
           
fc_test2    = [250 1000 2000 8000];
fmod_test2  = [20 30 50 70 100 150];
filenames2 = [];
for i = 1:length(fc_test2)
    for j = 1:length(fmod_test2)
        filenames2{i,j} = ['rough_test_fc_' Num2str(fc_test2(i),4) '_AM_m_100_fmod_' Num2str(fmod_test2(j),3) 'Hz'];
    end
end

fmod_test3  = [20 70 200];
filenames3  = {'track_39_t01', 'track_39_t04', 'track_39_t07'; ...
               'track_39_t02', 'track_39_t05', 'track_39_t08'; ...
               'track_39_t03', 'track_39_t06', 'track_39_t09'};
titles3     = {'AM BBN', 'AM Tone', 'FM Tone'};
    
fmod_test4  = [20 70 200];
SPL_test4   = [40:10:80];
filenames4  = {'track_39_t02',60; ...
               'track_39_t05',70; ...
               'track_39_t08',70};
titles4     = {'AM BBN', 'AM Tone', 'FM Tone'};

fmod_test5  = 70;
fdev_test5  = [30 60 120 150 400 600 800 1000];
fc_test5    = 1500;
filenames5  = [];
for j = 1:length(fdev_test5)
    filenames5{j} = ['rough_test_fc_' Num2str(fc_test5,4) '_FM_dev_' Num2str(fdev_test5(j),3) '_fmod_' Num2str(fmod_test5,3) 'Hz_60_dBSPL'];
end
        
fmod_test6  = [0 10 20 40 50 60 70 80 100 120 150 200 500];
fdev_test6  = 800;
fc_test6    = 1600;
filenames6  = [];
for j = 1:length(fmod_test6)
    filenames6{j} = ['rough_test_fc_' Num2str(fc_test6,4) '_FM_dev_' Num2str(fdev_test6,3) '_fmod_' Num2str(fmod_test6(j),3) 'Hz_60_dBSPL'];
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 1
    outs.filenames = [filenames0; filenames1(:,1); filenames2(:); filenames3(:); filenames4(:,1); filenames5(:); filenames6(:)];
    outs.dir_where = {pathaudio; pathaudio_src};
    return;
end

ExpNo = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp0
    
    filename = [pathaudio filenames0{1} '.wav'];
    listed_files{count_files} = sprintf('%.0f - %s',ExpNo,filename);
    count_files = count_files + 1;

    [x Fs] = Wavread(filename);

    starti = 1;
    insig = x(starti:starti + N-1);

    out = Roughness_offline_debug(insig,Fs,N, optsDebug); %No padding needed for off-line version
    R0(1) = out{1};

    disp(sprintf('Exp 0: R=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames0{1})));

end
ExpNo = ExpNo + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp1 
    % Fastl2007, Fig 11.1 (it considers 60 dB SPL)
    % Daniel1997, Fig. 5 (it considers 70 dB SPL)
    
    for k = 1:length(filenames1)
    
        filename = [pathaudio_src filenames1{k,1} '.wav'];
        listed_files{count_files} = sprintf('%.0f - %s',ExpNo,filename);
        count_files = count_files + 1;
    
        mx(k) = filenames1{k,2};
        
        [x Fs] = Wavread(filename);

        starti = 1;
        insig_tmp = buffer(x,N,0);
        dBFS(k) = rmsdb(x);
        
        for j = 1:3
            insig = insig_tmp(:,j);
            out = Roughness_offline_debug(insig,Fs,N, optsDebug); % insig at dB SPL
            R(j,k)  = out{1};
            dBSPL(k) = out{3};
            disp(sprintf('R=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames1{k})))
        end
        
    end

    figure;
    % plot(mx,mean(R)), grid on
    errorbar(mx*100,mean(R),std(R)), grid on
    hold on
    Rteo = 1.36*mx.^1.6;
    plot(mx*100,Rteo,'r--');
    xlabel('modulation index m [%]')
    ylabel('Roughness [asper]')
    h(end+1) = gcf;
    title('1 kHz AM tone, fmod = 70 Hz, 70 dB SPL')
end
ExpNo = ExpNo + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp2
    % Fastl2007, Fig 11.2 (it considers 60 dB SPL)
    % Daniel1997, Fig. 3 (it considers 60 dB SPL)
    
    for i = 1:length(fc_test2)
        for j = 1:length(fmod_test2)
    
            txtLegend{i} = ['f_c=' num2str(fc_test2(i)) '[Hz]'];
            filename = [pathaudio filenames2{i,j} '.wav'];
            listed_files{count_files} = sprintf('%.0f - %s',ExpNo,filename);
            count_files = count_files + 1;
            
            [x Fs] = Wavread(filename);
            
            dBFS2(i) = rmsdb(x);
            
            starti = 1;
            insig = x(starti:starti + N-1);
            t = ( 1:length(insig) )/Fs;     
            out = Roughness_offline_debug(insig,Fs,N, optsDebug); %No padding needed for off-line version
            
            R2(i,j)  = out{1};
            dBSPL2(i) = out{3}; % dBFS2(i) + 90;
            disp(sprintf('Exp 2: R=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames2)))

        end
    end
    
    figure;
    plot(fmod_test2,R2), grid on
    legend( txtLegend )
    xlabel('modulation frequency [Hz]')
    ylabel('Roughness [asper]')
    h(end+1) = gcf;
    
end
ExpNo = ExpNo + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp3
    % Fastl2007, Fig 11.3 (it considers XX dB SPL)
    figure;
    for j = 1:3
    
        for k = 1:3
            filename = [pathaudio_src filenames3{k,j} '.wav'];
            listed_files{count_files} = sprintf('%.0f - %s',ExpNo,filename);
            count_files = count_files + 1;
            
            [x Fs] = Wavread(filename);

            starti = 1;
            insig = x(starti:starti + N-1);
            
            if k == 2
                disp('');
                warning('in debug...');
            end
            
            out = Roughness_offline_debug(insig,Fs,N, optsDebug); %No padding needed for off-line version
            R3(j,k) = out{1};
        
            disp(sprintf('Exp 3: R=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames3{k})))
        end
        
        subplot(1,3,j)
        plot(fmod_test3,R3(j,:),'o--'), grid on
        if j == 2
            xlabel('modulation frequency [Hz]')
        end
        ylabel('Roughness [asper]')
        title( titles3{j} )
                
    end
    h(end+1) = gcf;
    
end
ExpNo = ExpNo + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp4
    
    % Fastl2007, Fig 11.4: modulation depth = 40 dB
    
    figure;
    for k = 1:3
        
        d_dB = SPL_test4 - filenames4{k,2};
        filename = [pathaudio_src filenames4{k,1} '.wav'];
        listed_files{count_files} = sprintf('%.0f - %s',ExpNo,filename);
        count_files = count_files + 1;
        
        [x Fs] = Wavread(filename);
        starti = 1;
            
        for j = 1:length(d_dB)
            
            insig = From_dB(d_dB(j))*x(starti:starti + N-1);
           
            out = Roughness_offline_debug(insig,Fs,N, optsDebug); %No padding needed for off-line version
            R4(k,j) = out{1};
        
            disp(sprintf('Exp 4: R=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames4{k})))
        end
        
        subplot(1,3,k)
        plot(SPL_test4,R4(k,:)), grid on
        if k == 2
            xlabel('Sound Pressure Level [dB]')
        end
        ylabel('Roughness [asper]')
        title( titles4{k} )
               
    end
    h(end+1) = gcf;
    
end
ExpNo = ExpNo + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp5
    % Fastl2007, Fig.11.5
    
    for j = 1:length(fdev_test5)
    
        filename = [pathaudio filenames5{j} '.wav'];
        listed_files{count_files} = sprintf('%.0f - %s',ExpNo,filename);
        count_files = count_files + 1;
        
        [x Fs] = Wavread(filename);

        starti = 1;
        insig = x(starti:starti + N-1);
        out = Roughness_offline_debug(insig,Fs,N, optsDebug); %No padding needed for off-line version

        R5(j)  = out{1};

        disp(sprintf('Exp 2: R=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames5{j})))

    end
    
    figure;
    plot(fdev_test5,R5), grid on
    xlabel('frequency deviation [Hz]')
    ylabel('Roughness [asper]')
    h(end+1) = gcf;
    title(sprintf('FM tone, fc=1500 [Hz], fmod=%.0f [Hz], 60 dB SPL',fmod_test5))
    
end
ExpNo = ExpNo + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoExp6
    % Daniel1997, Fig.9
    
    for j = 1:length(fmod_test6)
    
        filename = [pathaudio_D filenames6{j} '.wav'];
        listed_files{count_files} = sprintf('%.0f - %s',ExpNo,filename);
        count_files = count_files + 1;
        
        [x Fs] = Wavread(filename);

        starti = 1;
        insig = x(starti:starti + N-1);
        out = Roughness_offline_debug(insig,Fs,N, optsDebug); %No padding needed for off-line version

        R6(j)  = out{1};

        disp(sprintf('Exp 6: R=%.3f [aspers]\t test signal: %s\n',out{1},name2figname(filenames6{j})))

    end
    
    figure;
    plot(fmod_test6,R6), grid on
    xlabel('modulation frequency [Hz]')
    ylabel('Roughness [asper]')
    h(end+1) = gcf;
    title(sprintf('FM tone, fc=%.0f [Hz], dev.freq=%.0f, 60 dB SPL',fc_test6,fdev_test6))
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Save_all_figures(h,pathfigures);

for k = 1:length(listed_files)
    disp(listed_files{k})
end

if bDiary
	diary off
end

end

