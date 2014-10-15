function [info output] = experiment_report_20140228_Physical_validation_LIST(info)
% function [info output] = experiment_report_20140228_Physical_validation_LIST(info)
%
%     info.bAssess	1	assess physical measures
%     info.bAnalyse
%     info.root_dir	'~/Documenten/fda_eval_LIST/'
%     info.F0reference	'praat'
%     info.F0max	400
%     info.Fs		15659.375
%     info.bCleanSpeech 	0
%     info.typenoise    'SSN'
%
% Stand alone example:
%       options.typenoise = 'white';
%       experiment_report_20140228_Physical_validation_LIST(options);
% 
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2014
% Created on    : 28/02/2014
% Last update on: 18/09/2014 % Update this date manually
% Last use on   : 09/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;

if nargin == 0
    try
        localpaths = Get_paths([],1); % addpath('~/Documenten/MATLAB/MATLAB_svn/Utility/')
    catch
        warning('KU Leuven set-up (Get_paths) not found');
    end
    info = [];
end

info = Ensure_field(info, 'bAssess', 1); % To assess physical measures
info = Ensure_field(info, 'bAnalyse', 1);
info = Ensure_field(info, 'typenoise','SSN'); % 'white';
info = Ensure_field(info,'bSave', 1); % F0mod_validation_PB_database
info.tstep       = 5e-4;

if nargout > 1
    bPlot = 0;
else
    bPlot = 1;
end

close all
if isunix
    info  	= Ensure_field(info,'root_dir'      ,'~/Documenten/fda_eval_LIST/');
else
    try 
        tmp = Get_TUe_subpaths('db_speechmaterials'); % DELL computer
        
        info.wavinfo    = tmp.allfiles_LISTf;
        info.root_dir   = tmp.fda_eval_LISTf;
    catch
        root_dir = uigetdir;
        info.root_dir = [root_dir delim];
    end
end

if ~isfield(info,'bCleanSpeech')
    
    if strcmp(info.F0reference, 'praat')
        % info.bClearSpeech   = input('Press 1 if you want to use clean speech for Praat-F0 estimation. Press 0 if you want to use CP810 simulated speech: ');
        info.bCleanSpeech   = 0; % if 0 then emulated CP810 wav files will be used
    end
    
end

if info.bCleanSpeech == 0
    info = Ensure_field(info,'output',[info.root_dir 'output-CP810-' info.typenoise delim]);
else
    info = Ensure_field(info,'output',[info.root_dir 'output-' info.typenoise delim]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Diary(mfilename, bDiary, info.output);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.F0reference = 'praat';
info        = Ensure_field(info,'F0max'         , 400);
info        = Ensure_field(info, 'Fs'           , round(15659.375));
info.speaker = 'wdz';

p.CFG.Fs        = info.Fs;
p.F0mod.F0max   = info.F0max;

if strcmp(info.F0reference,'praat')
    
    if info.bCleanSpeech == 0
        
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Using audio files from CP810 simulation... (bCleanSpeech = 0)')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        info.praat = [info.root_dir 'praat-CP810' delim];   % /praat-CP810/ 

    else
        
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('Using original audio files... (bCleanSpeech = 1)')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        info = Ensure_field(info, 'praat', [info.root_dir 'praat'       delim]);
        
    end
    
    if info.bAssess
        bCreated = Mkdir(info.praat);
    end
end

if info.bAssess
    % Tested on 20/09/2014
    bCreated = Mkdir(info.output);

    if bCreated == 0
        warning(['Output directory non empty, please rename folder ' info.output ' and re-run this script'])

        disp('  ')
        disp('Press any button to continue');
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        % pause()
    end
    
    info = F0mod_validation_LISTf_database(info, p);
   
end

if strcmp(info.speaker,'sb') | strcmp(info.speaker,'wdz')
    idx_to_do = 1;
elseif strcmp(info.speaker,'rl')
    idx_to_do = 2;
end

Colors = [0 0 0; 1 1 1];

if info.bAnalyse
    
    speakers = {'wdz','rl'};
    gender   = {' (female)', ' (male)'};
    
    info.nSentences = 350;
    N = info.nSentences;
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp(sprintf('%s.m: %.0f files to be analysed',mfilename,N))
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    pause(2)
    
    for i = idx_to_do
        h = [];
        
        if info.bCleanSpeech == 0
            info = Ensure_field(info,'t_silence',0.5);
            t_silence = info.t_silence; 
            fprintf('%s.m: Compensating unvoiced time added in Simulink (150ms begin and 350ms at the end) in a total of %.2f [s]\n',mfilename,t_silence);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            pause(2)
        else
            t_silence = 0;
        end
        
        database_info = import_physical_measure([info.output speakers{i} '-database-time-info.txt'  ], 2, N+1); 
        t_voiced = database_info(:,3);
        t_unvoiced = database_info(:,4)-t_silence;
        t_total = database_info(:,2)-t_silence;
        vErrSim     = import_physical_measure([info.output speakers{i} '-vErrSim.txt'  ], 2, N+1);
        uvErrSim    = import_physical_measure([info.output speakers{i} '-uvErrSim.txt' ], 2, N+1);
        gErrSim     = import_physical_measure([info.output speakers{i} '-gErrSim.txt'  ], 2, N+1);
        f0DevErrSim = import_physical_measure([info.output speakers{i} '-f0ErrSim.txt' ], 2, N+1);
        
        numSNR = size(vErrSim,2)-1;
        [m1  s1 ] = prepare_barplot(vErrSim(:,2:end));
        [m11 s11] = prepare_barplot(vErrSim(:,2:end).*repmat(t_voiced./t_total,1,numSNR));
        
        stPlot.yLim         = [0 100];
        stPlot.YLabel       = 'vErr [%]';
        m1label = stPlot.YLabel;
        if bPlot
            figure
            h(end+1) = Plot_measure(m1,s1,[speakers{i} gender{i}], stPlot, Colors(i,:));
        end
        stPlot = [];
        
        [m2 s2] = prepare_barplot(uvErrSim(:,2:end));
        [m12 s12] = prepare_barplot(uvErrSim(:,2:end).*repmat(t_unvoiced./t_total,1,numSNR));
        
        stPlot.yLim         = [0 25];
        stPlot.YLabel       = 'uvErr [%]';
        m2label = stPlot.YLabel;
        if bPlot
            figure
            h(end+1) = Plot_measure(m2,s2,[speakers{i} gender{i}], stPlot, Colors(i,:));
        end
        stPlot = [];
        
        [m3 s3] = prepare_barplot( gErrSim(:,2:end));
        [m13 s13] = prepare_barplot(gErrSim(:,2:end).*repmat(t_voiced./t_total,1,numSNR));
        
        [m99 s99] = prepare_barplot(    vErrSim(:,2:end) .*repmat(t_voiced./t_total,1,numSNR) + ...
                                        uvErrSim(:,2:end).*repmat(t_unvoiced./t_total,1,numSNR) + ...
                                        gErrSim(:,2:end) .*repmat(t_voiced./t_total,1,numSNR));
        stPlot.yLim         = [0 25];
        stPlot.YLabel       = 'gErr [%]';
        m3label = stPlot.YLabel;
        if bPlot
            figure
            h(end+1) = Plot_measure(m3,s3,[speakers{i} gender{i}], stPlot, Colors(i,:));
        end
        stPlot = [];
        
        [m4 s4] = prepare_barplot(f0DevErrSim(:,2:end));
        figure
        stPlot.yLim         = [0 25];
        stPlot.YLabel       = 'f0Dev [\Delta Hz]';
        m4label = stPlot.YLabel;
        if bPlot
            h(end+1) = Plot_measure(m4,s4,[speakers{i} gender{i}], stPlot, Colors(i,:));
        end
        stPlot = [];
        
        output.m1 = m1;
        output.m1label = m1label;
        output.s1 = s1;
        output.m2 = m2;
        output.m2label = m2label;
        output.s2 = s2;
        output.m3 = m3;
        output.m3label = m3label;
        output.s3 = s3;
        output.m4 = m4;
        output.m4label = m4label;
        output.s4 = s4;
        
        output.m11 = m11;
        output.s11 = s11;
        output.m12 = m12;
        output.s12 = s12;
        output.m13 = m13;
        output.s13 = s13;
        output.m99 = m99;
        output.s99 = s99;
        
        output.t_voiced = sum(t_voiced);
        output.t_unvoiced = sum(t_unvoiced);
        output.t_total = sum(t_total);
        
        output.vErrSim = vErrSim;
        output.uvErrSim = uvErrSim;
        output.gErrSim = gErrSim;
        
        output.snr = [99 20 10 5 0 -5];

        
        if info.bSave & bPlot
            for k = 1:length(h)
                try
                    Saveas(h(k),[info.figures_folder    speakers{i} '-' num2str(k)]);
                catch
                    Saveas(h(k),[speakers{i} '-' num2str(k)]);
                end
            end
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDiary
	diary off
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('EOF')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m s] = prepare_barplot(x, y, z)
% function [m s] = prepare_barplot(x, y, z)
%
% m - mean
% s - std
m = [];
s = [];

if exist('x','var')
    [m1    ,s1    ]     = Get_mean(x);
%     [m1_tot,s1_tot]     = Get_mean(m1');
    m = [m m1'];
    s = [s s1'];
end
if exist('y','var')
    [m2    ,s2    ]     = Get_mean(y);
%     [m2_tot,s2_tot]     = Get_mean(m2');
    m = [m m2'];
    s = [s s2'];
end

if exist('z','var')
    [m3    ,s3    ]     = Get_mean(z);
%     [m3_tot,s3_tot]     = Get_mean(m3');
    m = [m m3'];
    s = [s s3'];
end

function h = Plot_measure(m,s,title_str,stPlot, Colors)

if nargin < 4
    stPlot = [];
end

nTicks = 10;
stPlot.figPos       = [0 0 1024 300];
stPlot.xTickLabel   = {'Q','20','10','5','0','-5'};
stPlot.xTick        = 1:length(stPlot.xTickLabel);
stPlot.xLim         = [ 0   max(stPlot.xTick)+1];
stPlot              = Ensure_field(stPlot, 'yLim', [-2 20]);
stepTicks           = ceil( (stPlot.yLim(2)-stPlot.yLim(1))/nTicks );

stPlot.yTick        = [stPlot.yLim(1)+stepTicks:stepTicks:stPlot.yLim(2)-stepTicks];
stPlot.Title        = title_str;
stPlot.XLabel       = 'SNR (dB)';
stPlot              = Ensure_field(stPlot, 'YLabel','%');

if size(m,2) == 1
    
    if ~exist('Colors','var')
        Colors              = [1 1 1];  
    end
    m = [m 0*m]; % trick to plot bars series with x-axis in steps of 1
    s = [s 0*s];
else % then size == 2
    stPlot.SeriesLabel  = {'NMT','Sim'};
    if ~exist('Colors','var')
        Colors              = [1 1 1; 0 0 0];  
    end
end

if isfield(stPlot, 'SeriesLabel')
    barweb(m,s,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors,[],stPlot.SeriesLabel);
else
    barweb(m,s,[],[],stPlot.Title,stPlot.XLabel,stPlot.YLabel,Colors);
end

grid on, hold on
plot(stPlot.xLim,[50 50],'r--')

h = gcf;
set(h,'Position', stPlot.figPos);

Handle = gca;
set(Handle,'YTick',stPlot.yTick)
set(Handle,'XTick',stPlot.xTick)
set(Handle,'XLim',stPlot.xLim)
set(Handle,'YLim',stPlot.yLim)
set(Handle,'XTickLabel',stPlot.xTickLabel)

set(h, 'PaperPositionMode','auto')

