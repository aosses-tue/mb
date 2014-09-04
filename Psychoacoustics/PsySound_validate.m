function h = PsySound_validate(options)
% function h = PsySound_validate(options)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       options.bDoFluct = 1;
%       PsySound_validate(options);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 17/08/2014
% Last update on: 17/08/2014 % Update this date manually
% Last use on   : 17/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = [];

if nargin == 0
    options = [];
end

if ~isunix
    dir_where_ref = 'D:\MATLAB\Output\tmp-cal\';
else
    dir_where_ref = Get_TUe_paths('outputs'); % tested on 17/08/2014
end

close all

bDiary = 0;
Diary(mfilename,bDiary,dir_where_ref);

options.bDoLoud = 0;
options.bDoSharp = 0;
options.bDoRough = 0;
options = Ensure_field(options,'bDoFluct',0);
bSave = 1;

if bSave == 1
    Mkdir([dir_where_ref 'Figures' delim]);
end

t1 = 0;
t2 = 4+t1;

bGenerateTestSounds = 0;
% use: Generate_reference_sounds to generate reference sounds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cal file: white noise:
options.calfile = [dir_where_ref 'track_03.wav'];
options.callevel = 60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref_loud  = [dir_where_ref 'ref_loud.wav'];
ref_sharp = [dir_where_ref 'ref_sharp.wav'];
ref_fluct = [dir_where_ref 'ref_fluct.wav'];
ref_rough = [dir_where_ref 'ref_rough.wav'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Loudness

if options.bDoLoud % validated on 18/08/2014
    
    hLoud = [];
    
    if bGenerateTestSounds
        [x fs] = Wavread(ref_loud); % 40 dB SPL
        for i = 10:10:30 % 40 to 70
            y = From_dB(i)*x;
            fprintf('RMS value of test signal = %.2f [dB SPL]. (0 dBFS = 90 dB SPL)\n',rmsdb(y)+90);

            Wavwrite(y, fs, [dir_where_ref 'test_loud_' num2str(40+i) 'dB']);
        end
    end

    option.bExtension = 0;
    loud_files =Get_filenames(dir_where_ref,'*_loud*.wav',option);

    for i = 1:length(loud_files)
        options.nAnalyser = 12; % DLM Chalupper
        [tmp_h tmp output] = PsySoundCL([dir_where_ref loud_files{i} '.wav'],options);
        close(tmp_h([1:2 4:6]))
        hLoud(end+1) = tmp_h(3);
        Mean_loud = mean(output.DataLoud( find(output.t > t1 & output.t < t2) ));
        Min_loud = prctile(output.DataSharp( find(output.t > t1 & output.t < t2) ), 5);
        Max_loud = prctile(output.DataSharp( find(output.t > t1 & output.t < t2) ),95);
        
        fprintf('Mean loudness (from t1 = %.3f to t2 = %.3f) = %.3f [sone]\n',t1,t2,Mean_loud);
        fprintf('Min and max values (percentiles 5 and 95 [%%] = %.3f and %.3f respectively)',Min_loud,Max_loud);
        
        if bSave == 1
            Print_date_on_figure;
            Saveas(hLoud(end),[dir_where_ref 'Figures' delim loud_files{i}]);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Sharpness

if options.bDoSharp
    hSharp = [];
    
    option.bExtension = 0;
    sharp_files =Get_filenames(dir_where_ref,'*_sharp*.wav',option);

    for i = 1:length(sharp_files)
        options.nAnalyser = 12; % DLM Chalupper
        [tmp_h tmp output] = PsySoundCL([dir_where_ref sharp_files{i} '.wav'],options);
        close(tmp_h([1:5]))
        hSharp(end+1) = tmp_h(6);
        Mean_sharp = Get_mean(output.DataSharp( find(output.t > t1 & output.t < t2) ));
        Min_sharp = prctile(output.DataSharp( find(output.t > t1 & output.t < t2) ),5);
        Max_sharp = prctile(output.DataSharp( find(output.t > t1 & output.t < t2) ),95);
        
        fprintf('Mean sharp (from t1 = %.3f to t2 = %.3f) = %.3f [acum]\n',t1,t2,Mean_sharp);
        fprintf('Min and max values (percentiles 5 and 95 [%%] = %.3f and %.3f respectively)\n',Min_sharp,Max_sharp);
        
        if bSave == 1
            Print_date_on_figure;
            Saveas(hSharp(end),[dir_where_ref 'Figures' delim sharp_files{i}]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Fluctuation strength
% 3.1 AM Tones
if options.bDoFluct

    dBFS = 100; 
    option.bExtension = 0;
    fc = 1000;
    fluct_files =Get_filenames(dir_where_ref,['*_fluct_fc_' num2str(fc) '*.wav'],option);
    % [f1 f2] = fileparts('D:\Databases\dir04-Psychoacoustics\Fastl-and-Zwicker-2007\03-Extracted-files\track_35_t05.wav')
    % fluct_files{1} = f2;
    % dir_where_ref = [f1 delim];
    % fmod = 4;
    lf = [];
    le_min = [];
    le_max = [];
    fmod = [0.5 1 2 4 8 16 32];
    for i=1:length(fluct_files)
        [x fs] = Wavread([dir_where_ref fluct_files{i} '.wav']);
        outstmp = Do_fluct(x,fs,dBFS);
        lf      = [lf; outstmp.lf];
        le_min  = [le_min; outstmp.le_min];
        le_max  = [le_max; outstmp.le_max];
        FS(i)   = ( 0.008*sum(outstmp.le_diff) )/( fmod(i)/4 + 4/fmod(i) );
        FS2(i)   = ( 0.008*sum(outstmp.le_max-outstmp.le_min) )/( fmod(i)/4 + 4/fmod(i) );
        FS2 - FS
    end
    
    if fmod(2)==1 & fmod(4)==4 & fmod(6)==16
        
        figure;
        subplot(1,2,1)
        plot(1:24, le_max([1 4],:),'o-','LineWidth',1), grid on, hold on
        plot(1:24, le_min([1 4],:),'LineWidth',3) 
        ylabel('L_E [dB]')
        xlabel('Critical-band rate [Bark]')
        ha = gca;
        legend('fmod = 1 Hz','fmod = 4 Hz')
        title(sprintf('AM Sine tones at %.0f [Hz]',fc))
        
        subplot(1,2,2)
        plot(1:24, le_max([6 4],:),'o-','LineWidth',1), grid on, hold on
        plot(1:24, le_min([6 4],:),'LineWidth',3) 
        ylabel('L_E [dB]')
        xlabel('Critical-band rate [Bark]')
        ha(end+1) = gca;
        legend('fmod = 16 Hz','fmod = 4 Hz')
        linkaxes(ha,'xy');
        xlim([5 20])
        h(end+1) = gcf;
        ylim([-10 90])
    end
    
    figure;
    semilogx(fmod,100*FS/max(FS),'bo-'), grid on
    xlim([0.5 25])
    xlabel('Modulation frequency [Hz]')
    ylabel('Relative fluctuation [%]')
    legend(sprintf('max(FS) = %.2f [vacil]',max(FS)))
    h(end+1) = gcf;
    title('FS for AM tones')
    % disp Fig 1: min and max Excitation patterns
    % disp Fig 2: Relative fluctuation strength
    ylim([0 110])
end

% 3.2 FM Tones
if options.bDoFluct

    option.bExtension = 0;
    fc = 1500;
    fluct_files =Get_filenames(dir_where_ref,['*_fluct_fc_' num2str(fc) '*.wav'],option);
    lf = [];
    le_min = [];
    le_max = [];
    fmod = [0.5 1 2 4 8 16 32];
    for i=1:length(fluct_files)
        [x fs] = Wavread([dir_where_ref fluct_files{i} '.wav']);
        x = From_dB(-10)*x;
        fprintf('RMS to be considered: %.2f dB SPL (0 dBFS = 110 dB SPL)\n',rmsdb(x)+110);
        outstmp = Do_fluct(x,fs);
        lf      = [lf; outstmp.lf];
        le_min  = [le_min; outstmp.le_min];
        le_max  = [le_max; outstmp.le_max];
        FS(i)   = ( 0.008*sum(outstmp.le_diff) )/( fmod(i)/4 + 4/fmod(i) );
    end
    
    if fmod(2)==1 & fmod(4)==4 & fmod(6)==16
        
        figure;
        subplot(1,2,1)
        plot(1:24, le_max([1 4],:),'o-','LineWidth',1), grid on, hold on
        plot(1:24, le_min([1 4],:),'LineWidth',3) 
        ylabel('L_E [dB]')
        xlabel('Critical-band rate [Bark]')
        ha = gca;
        legend('fmod = 1 Hz','fmod = 4 Hz')
        title(sprintf('FM Sine tones at %.0f [Hz]',fc))
        
        subplot(1,2,2)
        plot(1:24, le_max([6 4],:),'o-','LineWidth',1), grid on, hold on
        plot(1:24, le_min([6 4],:),'LineWidth',3) 
        ylabel('L_E [dB]')
        xlabel('Critical-band rate [Bark]')
        ha(end+1) = gca;
        legend('fmod = 16 Hz','fmod = 4 Hz')
        linkaxes(ha,'xy');
        xlim([5 20])
        h(end+1) = gcf;
        ylim([0 90])
    end
    
    figure;
    semilogx(fmod,100*FS/max(FS),'bo-'), grid on
    xlim([0.5 25])
    xlabel('Modulation frequency [Hz]')
    ylabel('Relative fluctuation [%]')
    legend(sprintf('max(FS) = %.2f [vacil]',max(FS)))
    h(end+1) = gcf;
    % disp Fig 3: min and max Excitation patterns
    % disp Fig 4: Relative fluctuation strength
    title('FS for FM tones')
    ylim([0 110])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Roughness

if options.bDoRough
    
    hRough = [];
    
    option.bExtension = 0;
    rough_files =Get_filenames(dir_where_ref,'*_rough*.wav',option);
    cont = 1;
    cont125 = 1;
    cont500 = 1;
    
    for i = 1:length(rough_files)
        options.nAnalyser = 15; % 
        [tmp_h tmp output] = PsySoundCL([dir_where_ref rough_files{i} '.wav'],options);
        
        close(tmp_h(1))
        hRough(end+1) = tmp_h(2);
        Mean_rough = Get_mean(output.DataRough( find(output.t > t1 & output.t < t2) ));
        Min_rough = prctile(output.DataRough( find(output.t > t1 & output.t < t2) ),5);
        Max_rough = prctile(output.DataRough( find(output.t > t1 & output.t < t2) ),95);
        
        fprintf('Mean roughness (from t1 = %.3f to t2 = %.3f) = %.3f [asper]\n',t1,t2,Mean_rough);
        fprintf('Min and max values (percentiles 5 and 95 [%%] = %.3f and %.3f respectively)\n',Min_rough,Max_rough);
        
        if Match_in_str('fc_1000_AM_m_100',rough_files{i}) % then is not the reference
            rough_value(cont) = Mean_rough;
            cont = cont + 1;
        end
       
        if Match_in_str('fc_0125_AM_m_100',rough_files{i}) % then is not the reference
            rough_value125(cont125) = Mean_rough;
            cont125 = cont125 + 1;
        end
        
        if Match_in_str('fc_0500_AM_m_100',rough_files{i}) % then is not the reference
            rough_value500(cont500) = Mean_rough;
            cont500 = cont500 + 1;
        end
        
        if bSave == 1
            Print_date_on_figure;
            Saveas(hRough(end),[dir_where_ref 'Figures' delim rough_files{i}]);
        end
    end
    
    output.f = 10:20:170;
    output.rough_value(1,:) = rough_value;
    output.rough_value(2,:) = rough_value125;
    output.rough_value(3,:) = rough_value500;
    
    figure;
    plot(output.f,output.rough_value(1,:),'o-','LineWidth',2), hold on
    plot(output.f,output.rough_value(2,:),'ro--')
    plot(output.f,output.rough_value(3,:),'kx-')
    
    grid on;
    xlabel('Frequency [Hz]')
    ylabel('Roughness [asper]')
    legend('f_c = 1 kHz', 'f_c = 125 Hz', 'f_c = 500 Hz')
    hRough(end+1) = gcf;
    ylim([-0.18 1.4])
    if bSave == 1
        Print_date_on_figure;
        Saveas(hRough(end),[dir_where_ref 'Figures' delim 'roughness-vs-freq']);
    end
    
end
if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
