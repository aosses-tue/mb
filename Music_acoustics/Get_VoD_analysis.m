function [output h] = Get_VoD_analysis(filename1,filename2,options)
% function [output h] = Get_VoD_analysis(filename1,filename2,options)
%
% 1. Description:
%       Same than VoD_comparisons, but here analysis is done over aligned/
%       truncated data.
% 
%       Make sure you have run previously the script VoD_run.m (without 
%       output parameters) in order to generate up-to-date calibrated 
%       wav-files
% 
% 2. Additional info:
%       Tested cross-platform: No
%       near-field: YES, OK
%       far-field : NO,  OK 
%
% 3. Stand-alone example:
%       Get_VoD_analysis;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/08/2014
% Last update on: 01/10/2014 % Update this date manually
% Last use on   : 01/10/2014 % Update this date manually
% 
% Original file name: VoD_comparisons2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    close all
    filename1 = [Get_TUe_paths('outputs') 'ref_loud.wav'];
    filename2 = [Get_TUe_paths('outputs') 'ref_rough.wav'];
    options.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav'];
    options.callevel = 60; % rms 90 dB SPL = 0 dBFS 
    
end
    
options     = Ensure_field(options,'label1',' meas');
options     = Ensure_field(options,'label2','model');

options     = Ensure_field(options,'bSave',0);
options     = Ensure_field(options,'bPlot',1);
options     = Ensure_field(options,'label','');
options     = Ensure_field(options,'SPLrange',[10 70]);
options     = Ensure_field(options,'frange',[50 5000]);

options     = Ensure_field(options,'bDoAnalyser08',0);
options     = Ensure_field(options,'bDoAnalyser10',0);
options     = Ensure_field(options,'bDoAnalyser11',0);
options     = Ensure_field(options,'bDoAnalyser12',1);
options     = Ensure_field(options,'bDoAnalyser15',0);

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

options.bCosineRamp = 0; % Cos ramp not yet applied for Loudness calculations
options.bCosineRampOnset = 0; %ms

options.bPlot = 0;

if options.bDoAnalyser08 == 1
    options.nAnalyser = 8;
    [tmp_h tmp_ha out_1_08]   = PsySoundCL(filename1,options);
    [tmp_h tmp_ha out_2_08]   = PsySoundCL(filename2,options);

    % % Example to plot SLM outputs
    %
    %     figure; 
    %     plot(   out_1_08.t,out_1_08.Data_dBAS, ...
    %             out_1_08.t,out_1_08.Data_dBAF, ...
    %             out_1_08.t,out_1_08.Data_dBZS, ...
    %             out_1_08.t,out_1_08.Data_dBZF); 
    %     legend('dB(A) S','dB(A) F', 'dB(Z) S', 'dB(Z) F')
    %     xlabel('Time (s)')
    %     ylabel('level (dB)')
    %     grid on
    %     title('SPL for audio file 1');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    options.nAnalyser = 15;
    [out_1_15 tmp_h tmp_ha] = PsySoundCL(filename1,options); % 2 figures
    [out_2_15 tmp_h tmp_ha] = PsySoundCL(filename2,options); % 2 figures
end

% % nAnalyser = 12
% 1. Average main loudness (Bark)       4. Main loudness (s)
% 2. Average specific loudness (Bark)   5. Specific loudness (s)
% 3. Loudness (s)                       6. Sharpness (s)
%
% % nAnalyser = 15
% 7. Specific roughness [Bark]
% 8. Roughness [s]

param = [];

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
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_12,out_2_12,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    
    if isfield(options,'ylim_loudness')
        ylim(options.ylim_loudness);
    end
    
    if isfield(options,'ylim_bExtend')
        
        if options.ylim_bExtend == 1
            y_old = get(gca,'YLim');
            x_old = get(gca,'XLim');
            ylim_extend(gca,1.25);
            plot([x_old],[y_old(2) y_old(2)],'k'); % horizontal line
        end
        
    end
    
    % Specific loudness
    options.nAnalyser   = 12;
    param{end+1}        = 'specific-loudness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_12,out_2_12,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    
    if isfield(options,'ylim_specific_loudness')
        ylim(options.ylim_specific_loudness);
    end
    
    if isfield(options,'ylim_bExtend')
        
        if options.ylim_bExtend == 1
            y_old = get(gca,'YLim');
            x_old = get(gca,'XLim');
            ylim_extend(gca,1.25);
            plot([x_old],[y_old(2) y_old(2)],'k'); % horizontal line
        end
        
    end
    
    % Sharpness
    options.nAnalyser   = 12;
    param{end+1}        = 'sharpness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_12,out_2_12,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.bDoAnalyser15 == 1
    options.nAnalyser = 15;
    param{end+1}        = 'roughness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_15,out_2_15,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));

    options.nAnalyser = 15;
    param{end+1} = 'specific-roughness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_15,out_2_15,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    
    if isfield(options,'ylim_specific_roughness')
        ylim(options.ylim_specific_roughness);
    end
        
    if isfield(options,'ylim_bExtend')
        
        if options.ylim_bExtend == 1
            y_old = get(gca,'YLim');
            x_old = get(gca,'XLim');
            ylim_extend(gca,1.25);
            plot([x_old],[y_old(2) y_old(2)],'k'); % horizontal line
        end
        
    end
    
end

try
%%%%
% Percentiles for Loudness
one_period_s        = diff(options.tanalysis) / (options.N_periods2analyse+1);
one_period_in_samples = ceil( one_period_s/min(diff(out_1_12.t)) );
% N_periods           = floor(options.time2save / one_period_s); % To analyse whole file
N_periods = options.N_periods2analyse+1;

pLmeas  = Get_percentiles_per_period(out_1_12.DataLoud,one_period_in_samples);
pLmodel = Get_percentiles_per_period(out_2_12.DataLoud,one_period_in_samples);

output.pL_in1   = pLmeas;
output.pL_in2   = pLmodel;
output.L_in1    = out_1_12.DataLoud;
output.L_in2    = out_2_12.DataLoud;
output.Lt       = out_1_12.t;

%%%%
% Percentiles for Roughness

one_period_s        = diff(options.tanalysis) / (options.N_periods2analyse+1);
one_period_in_samples = ceil( one_period_s/min(diff(out_1_15.t)) );
% N_periods           = floor(options.time2save / one_period_s);


pRmeas  = Get_percentiles_per_period(out_1_15.DataRough,one_period_in_samples);
pRmodel = Get_percentiles_per_period(out_2_15.DataRough,one_period_in_samples);

output.R_in1    = out_1_15.DataRough;
output.R_in2    = out_2_15.DataRough;
output.pR_in1   = pRmeas;
output.pR_in2   = pRmodel;
output.Rt       = out_1_15.t;

%%%%
% Percentiles for Sharpness
one_period_in_samples = ceil( one_period_s/min(diff(out_1_12.t)) );
yy1 = buffer(out_1_12.DataSharp,one_period_in_samples,0);
yy2 = buffer(out_2_12.DataSharp,one_period_in_samples,0);

y = yy1(:,1:N_periods);
p5 = percentile(y,5);
p50 = percentile(y,50);
p95 = percentile(y,95);

p5_mod = percentile(y,5);
p50_mod = percentile(y,50);
p95_mod = percentile(y,95);

catch
    warning('Percentile calculation not succeeded, maybe not every Analyser is enabled')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output.h = h; % figure handles

if options.bSave
    
    disp('Figures are going to be stored... Press ctr+C to cancel')
    pause(2)
    
    try
        paths.outputs = options.dst_folder_fig;
    catch
        paths.outputs   = Get_TUe_paths('outputs');
    end
    
    for i = 1:length(h)
        % options.format = 'emf';
        % Saveas(h(i),[options.dest_folder_fig 'psycho-fig-' param{i} options.label],options);
        options.format = 'epsc';
        Saveas(h(i),[options.dest_folder_fig 'psycho-fig-' param{i} options.label],options);
    end
    
else
    
    disp('Figures are NOT going to be stored... Set options.bSave to 1 and re-run the scripts in case you want to save the figures')
    pause(2)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end