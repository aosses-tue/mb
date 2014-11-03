function [output h] = Get_PsySound_analysis_3_files(filename1,filename2,filename3,options,stPlot)
% function [output h] = Get_PsySound_analysis_3_files(filename1,filename2,filename3,options,stPlot)
%
% 1. Description:
% 
% 2. Stand-alone example:
%       Get_PsySound_analysis;
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Original file : Get_VoD_analysis.m
% Created on    : 31/10/2014
% Last update on: 31/10/2014 % Update this date manually
% Last use on   : 31/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    close all
    filename1 = [Get_TUe_paths('outputs') 'ref_loud.wav'];
    filename2 = [Get_TUe_paths('outputs') 'ref_rough.wav'];
    options.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav'];
    options.callevel = 60; % rms 90 dB SPL = 0 dBFS     
end
    
if nargin < 4
    stPlot = [];
end

stPlot = Ensure_field(stPlot,'Title','PsySound Analysis');

options     = Ensure_field(options,'label1',' meas');
options     = Ensure_field(options,'label2','model');
options     = Ensure_field(options,'label3','meas*');

options     = Ensure_field(options,'bSave',0);
options     = Ensure_field(options,'bPlot',1);
options     = Ensure_field(options,'label','');
options     = Ensure_field(options,'SPLrange',[10 70]);
options     = Ensure_field(options,'frange',[50 5000]);
options     = Ensure_field(options,'frange_FFT',[50 5000]);

options     = Ensure_field(options,'bDoAnalyser01',0);
options     = Ensure_field(options,'bDoAnalyser08',0);
options     = Ensure_field(options,'bDoAnalyser10',0);
options     = Ensure_field(options,'bDoAnalyser11',0);
options     = Ensure_field(options,'bDoAnalyser12',1);
options     = Ensure_field(options,'bDoAnalyser15',0);

stPlot      = Ensure_field(stPlot,'bPlotHorLine',0);

if options.bSave == 1
    options = Ensure_field(options,'dest_folder_fig',Get_TUe_paths('outputs'));
end
    
h = []; % handles figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([mfilename '.m: analysing ' filename1 ' and ' filename2]);
 
ha1 = [];
ha2 = [];    
tmp_h = [];
tmp_h = [];

options = Ensure_field(options,'calfile',[Get_TUe_paths('db_calfiles') 'track_03.wav']);
options = Ensure_field(options,'callevel',70); % 'AMT' reference

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.bPlot = 0;

if options.bDoAnalyser01 == 1
    options.nAnalyser = 1;
    [out_1_01 tmp_h tmp_ha]   = PsySoundCL(filename1,options);
    [out_2_01 tmp_h tmp_ha]   = PsySoundCL(filename2,options);
end

if options.bDoAnalyser08 == 1
    options.nAnalyser = 8;
    [out_1_08 tmp_h tmp_ha]   = PsySoundCL(filename1,options);
    [out_2_08 tmp_h tmp_ha]   = PsySoundCL(filename2,options);
end

if options.bDoAnalyser10 == 1
    options.nAnalyser = 10;
    [out_1_10 tmp_h tmp_ha ]   = PsySoundCL(filename1,options);
    [out_2_10 tmp_h tmp_ha]   = PsySoundCL(filename2,options);
end

if options.bDoAnalyser11 == 1
    options.nAnalyser = 11;
    [out_1_11 tmp_h tmp_ha]   = PsySoundCL(filename1,options);
    [out_2_11 tmp_h tmp_ha]   = PsySoundCL(filename2,options);
end

if options.bDoAnalyser12 == 1
    options.nAnalyser = 12;
    [out_1_12 tmp_h tmp_ha]   = PsySoundCL(filename1,options); 
    [out_2_12 tmp_h tmp_ha]   = PsySoundCL(filename2,options); 
    [out_3_12 tmp_h tmp_ha]   = PsySoundCL(filename3,options); 
end

if options.bDoAnalyser15 == 1
    options.nAnalyser = 15;
    [out_1_15 tmp_h tmp_ha] = PsySoundCL(filename1,options); % 2 figures
    [out_2_15 tmp_h tmp_ha] = PsySoundCL(filename2,options); % 2 figures
end

param = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.bDoAnalyser01 == 1
    options.nAnalyser   = 1;
    param{end+1}        = 'FFT';
    [h(end+1) ha stats] = PsySoundCL_Figures(param{end},out_1_01,out_2_01,options,stPlot);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    
    strTitle = sprintf('Spectrum %s', stPlot.Title);
    title(strTitle);
    
    grid on
    xlim(options.frange_FFT);
    
    freqrange = options.frange_FFT(2)-options.frange_FFT(1);
    XTicks = round(options.frange_FFT(1):freqrange/4:options.frange_FFT(2));
    
    set(ha(end),'XTick',XTicks);
    set(ha(end),'YLim',options.SPLrange);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.bDoAnalyser10 == 1
    options.nAnalyser   = 10;
    param{end+1}        = 'specific-loudness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_10,out_2_10,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.bDoAnalyser11 == 1
    options.nAnalyser   = 11;
    param{end+1}        = 'one-third-OB';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_11,out_2_11,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    
    if length(stats.idx) == length(stats.t)
        legend( sprintf('tot = %.2f dB(Z)',dbsum( out_1_11.DataSpecOneThirdAvg )) ,...
                sprintf('tot = %.2f dB(Z)',dbsum( out_2_11.DataSpecOneThirdAvg )) );
    else
        legend( sprintf('tot = %.2f dB(Z)',dbsum( dbmean( out_1_11.DataSpecOneThird(stats.idx,:) )) ),...
                sprintf('tot = %.2f dB(Z)',dbsum( dbmean( out_2_11.DataSpecOneThird(stats.idx,:) )) ) );
    end
    
    if isfield(options,'SPLrange')
        set(gca,'YLim',options.SPLrange);
        
        plot(options.frange,[65 65],'k'); % horizontal line
    end
    
    if isfield(options,'frange')
        set(gca,'XLim',options.frange);
    end
    hold on
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.bDoAnalyser12 == 1
    
    % Loudness
    options.nAnalyser   = 12;
    param{end+1}        = 'loudness';
    options.color{3} = 'k-'; 
    [h(end+1) xx stats] = PsySoundCL_Figures_3(param{end},out_1_12,out_2_12,out_3_12,options,stPlot);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    
    strTitle = sprintf('Loudness %s', stPlot.Title);
    title(strTitle);
    
    if isfield(stPlot,'ylim_loudness')
        ylim(stPlot.ylim_loudness);
    end
    
    if isfield(stPlot,'ylim_bExtend')
        
        if stPlot.ylim_bExtend == 1
            y_old = get(gca,'YLim');
            x_old = get(gca,'XLim');
            ylim_extend(gca,1.25);
            if stPlot.bPlotHorLine == 1
                plot([x_old],[y_old(2) y_old(2)],'k'); % horizontal line
            end
        end
        
    end
    
    % Specific loudness
    options.nAnalyser   = 12;
    param{end+1}        = 'specific-loudness';
    [h(end+1) xx stats] = PsySoundCL_Figures_3(param{end},out_1_12,out_2_12,out_3_12,options,stPlot);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    
    strTitle = sprintf('Specific loudness %s', stPlot.Title);
    title(strTitle);
    
    if isfield(stPlot,'ylim_specific_loudness')
        ylim(stPlot.ylim_specific_loudness);
    end
    
    if isfield(stPlot,'ylim_bExtend')
        
        if stPlot.ylim_bExtend == 1
            y_old = get(gca,'YLim');
            x_old = get(gca,'XLim');
            ylim_extend(gca,1.25);
            if stPlot.bPlotHorLine == 1 
                plot([x_old],[y_old(2) y_old(2)],'k'); % horizontal line
            end
        end
        
    end
    
    % Sharpness
    options.nAnalyser   = 12;
    param{end+1}        = 'sharpness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_12,out_2_12,options,stPlot);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.bDoAnalyser15 == 1
    options.nAnalyser = 15;
    param{end+1}        = 'roughness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_15,out_2_15,options,stPlot);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));

    if isfield(options,'ylim_roughness')
        ylim(options.ylim_roughness);
    end
    
    options.nAnalyser = 15;
    param{end+1} = 'specific-roughness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_15,out_2_15,options,stPlot);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    
    if isfield(stPlot,'ylim_specific_roughness')
        ylim(stPlot.ylim_specific_roughness);
    end
        
    if isfield(stPlot,'ylim_bExtend')
        
        if stPlot.ylim_bExtend == 1
            y_old = get(gca,'YLim');
            x_old = get(gca,'XLim');
            ylim_extend(gca,1.25);
            plot([x_old],[y_old(2) y_old(2)],'k'); % horizontal line
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(h)
    hname{i} = ['psycho-fig-' param{i} options.label];
end
    
output.h = h; % figure handles
output.hname = hname;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end