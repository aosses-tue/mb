function h = Get_VoD_analysis(filename1,filename2,options)
% function h = Get_VoD_analysis(filename1,filename2,options)
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
% Last update on: 04/09/2014 % Update this date manually
% Last use on   : 04/09/2014 % Update this date manually
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
    
bDoChalupper    = 1;

% bHPF        = 1; % 1 = to load band-pass filtered audio files
options     = Ensure_field(options,'bSave',0);
options.bPlot  = 1;
options     = Ensure_field(options,'label','');

if options.bSave == 1
    options = Ensure_field(options,'dest_folder_fig',Get_TUe_paths('outputs'));
end
    
h = []; % handles figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([mfilename '.m: analysing ' filename1 ' and ' filename2]);
 
if bDoChalupper
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
    
    options.nAnalyser = 8;
    [tmp_h tmp_ha out_1_08]   = PsySoundCL(filename1,options);
    [tmp_h tmp_ha out_2_08]   = PsySoundCL(filename2,options);
%     [tmp_h tmp_ha out_1_10]   = PsySoundCL(filename1,options);
%     [tmp_h tmp_ha out_2_10]   = PsySoundCL(filename2,options);
    
    if options.bPlot
        figure; 
        plot(   out_1_08.t,out_1_08.Data_dBAS, ...
                out_1_08.t,out_1_08.Data_dBAF, ...
                out_1_08.t,out_1_08.Data_dBZS, ...
                out_1_08.t,out_1_08.Data_dBZF); 
        legend('dB(A) S','dB(A) F', 'dB(Z) S', 'dB(Z) F')
        xlabel('Time (s)')
        ylabel('level (dB)')
        grid on
        title('SPL for audio file 1');
    end
    options.nAnalyser = 10;
    [tmp_h tmp_ha out_1_10]   = PsySoundCL(filename1,options);
    [tmp_h tmp_ha out_2_10]   = PsySoundCL(filename2,options);
    
    options.nAnalyser = 11;
    [tmp_h tmp_ha out_1_11]   = PsySoundCL(filename1,options);
    [tmp_h tmp_ha out_2_11]   = PsySoundCL(filename2,options);
    
    options.nAnalyser = 12;
    [tmp_h tmp_ha out_1_12]   = PsySoundCL(filename1,options); % , times2zoom); % 6 figures
    [tmp_h tmp_ha out_2_12]   = PsySoundCL(filename2,options); % , times2zoom); % 6 figures
    
    options.nAnalyser = 15;
    [tmp_h tmp_ha out_1_15] = PsySoundCL(filename1,options); % 2 figures
    [tmp_h tmp_ha out_2_15] = PsySoundCL(filename2,options); % 2 figures
    
    % % nAnalyser = 12
    % 1. Average main loudness (Bark)       4. Main loudness (s)
    % 2. Average specific loudness (Bark)   5. Specific loudness (s)
    % 3. Loudness (s)                       6. Sharpness (s)
    %
    % % nAnalyser = 15
    % 7. Specific roughness [Bark]
    % 8. Roughness [s]
  
    param = [];
    
    options.nAnalyser   = 10;
    param{end+1}        = 'one-third-OB';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_10,out_2_10,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    legend('meas','model');
    
    options.nAnalyser   = 10;
    param{end+1}        = 'specific-loudness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_10,out_2_10,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    % legend('meas','model');
    legend( sprintf(' meas, tot = %.2f [sone]',out_1_10.stats.loud_tot), ... 
            sprintf('model, tot = %.2f [sone]',out_2_10.stats.loud_tot) );
    
    options.nAnalyser   = 11;
    param{end+1}        = 'one-third-OB';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_11,out_2_11,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    legend('meas','model');
    
    options.nAnalyser   = 12;
    param{end+1}        = 'loudness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_12,out_2_12,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    legend('meas','model');
    
    options.nAnalyser   = 12;
    param{end+1}        = 'sharpness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_12,out_2_12,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    legend('meas','model');
    
    options.nAnalyser = 15;
    param{end+1}        = 'roughness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_15,out_2_15,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    legend('meas','model');
    
    options.nAnalyser   = 12;
    param{end+1}        = 'specific-loudness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_12,out_2_12,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    legend('meas','model');
    
    options.nAnalyser = 15;
    param{end+1} = 'specific-roughness';
    [h(end+1) xx stats] = PsySoundCL_Figures(param{end},out_1_15,out_2_15,options);
    param{end} = sprintf('%s-analyser-%s',param{end},Num2str(options.nAnalyser));
    legend('meas','model');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        % Saveas(h(i),[options.dest_folder_fig 'psycho-fig' param{i} options.label],options);
        options.format = 'epsc';
        Saveas(h(i),[options.dest_folder_fig 'psycho-fig' param{i} options.label],options);
    end
    
else
    
    disp('Figures are NOT going to be stored... Set options.bSave to 1 and re-run the scripts in case you want to save the figures')
    pause(2)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end