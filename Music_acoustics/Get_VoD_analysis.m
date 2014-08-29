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
%       bDiary = 1; % to generate a log-file
%       VoD_comparisons2(bDiary);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 20/08/2014
% Last update on: 20/08/2014 % Update this date manually
% Last use on   : 20/08/2014 % Update this date manually
% 
% Original file name: VoD_comparisons2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDoChalupper    = 1;

bHPF        = 1; % 1 = to load band-pass filtered audio files
options.bSave  = 1;
options.bPlot  = 1;
options = Ensure_field(options,'label','');

if options.bSave == 1
    options = Ensure_field(options,'dest_folder_fig',Get_TUe_paths('outputs'));
end
    
h = []; % handles figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([mfilename '.m: analysing ' filename1 ' and ' filename2]);
 
% 4. Chalupper's model
%       - Last used on: 01/08/2014

if bDoChalupper
    ha1 = [];
    ha2 = [];    
    tmp_h = [];
    tmp_h = [];
    
    options.calfile = [Get_TUe_paths('db_calfiles') 'track_03.wav'];
    options.callevel = 70; % 'AMT' reference

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    options.bCosineRamp = 0; % Cos ramp not yet applied for Loudness calculations
    options.bCosineRampOnset = 0; %ms
    options.nAnalyser = 12;
    [tmp_h tmp_ha out_summary1]   = PsySoundCL(filename1,options); % , times2zoom); % 6 figures
    [tmp_h tmp_ha out_summary2]   = PsySoundCL(filename2,options); % , times2zoom); % 6 figures
    
    options.nAnalyser = 15;
    [tmp_h tmp_ha out_summary1_r] = PsySoundCL(filename1,options); % 2 figures
    [tmp_h tmp_ha out_summary2_r] = PsySoundCL(filename2,options); % 2 figures
    
    % % nAnalyser = 12
    % 1. Average main loudness (Bark)       4. Main loudness (s)
    % 2. Average specific loudness (Bark)   5. Specific loudness (s)
    % 3. Loudness (s)                       6. Sharpness (s)
    %
    % % nAnalyser = 15
    % 7. Specific roughness [Bark]
    % 8. Roughness [s]
    
    [h(end+1) xx stats] = PsySoundCL_Figures('loudness' ,out_summary1  ,out_summary2  ,options);
    legend('meas','model');
    
    [h(end+1) xx stats] = PsySoundCL_Figures('sharpness',out_summary1  ,out_summary2  ,options);
    % legend('meas','model');
    
    [h(end+1) xx stats] = PsySoundCL_Figures('roughness',out_summary1_r,out_summary2_r,options);
    % legend('meas','model');
    
    [h(end+1) xx stats] = PsySoundCL_Figures('specific-loudness' ,out_summary1  ,out_summary2  ,options);
    
    [h(end+1) xx stats] = PsySoundCL_Figures('specific-roughness',out_summary1_r,out_summary2_r,options);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if options.bSave
    
    disp('Figures are going to be stored...press any button to continue. Press ctr+C to cancel')
    pause()
    
    try
        paths.outputs = options.dst_folder_fig;
    catch
        paths.outputs   = Get_TUe_paths('outputs');
    end
    
    for j = 1:length(h)
        options.format = 'emf';
        Saveas(h(j),[options.dest_folder_fig 'psycho-fig' Num2str(j,2) options.label],options);
        options.format = 'epsc';
        Saveas(h(j),[options.dest_folder_fig 'psycho-fig' Num2str(j,2) options.label],options);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end