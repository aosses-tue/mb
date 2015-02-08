function r20150206_update
% function r20150206_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 06/02/2015
% Last update on: 06/02/2015 % Update this date manually
% Last use on   : 06/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bLearningCombFilter = 1;

if ~isunix
    dir ='D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-02-06-update\Hummer-audio-files\';
else
    
end

% Part 3
if bLearningCombFilter
    
    %% Information:
    c = 340; % in m/s
    K = 8192*4/2; % to plot FFT
    alfa = 1;
    fmaxPlot = 1200;
    
    % Initial position:
    S1   = [0 0  2.23];
    S1im = [0 0 -2.23];
    S2   = [0,.67, 2.23];
    S2im = [0,.67,-2.23];
    
    MicD = [1.58,0,1.68];
    
    distS1_S2       = norm(S2-S1);
    distS1_MicD     = norm(S1-MicD);
    distS1im_MicD   = norm(S1im-MicD);
    distS2_MicD     = norm(S2-MicD);
    distS2im_MicD   = norm(S2im-MicD);
    
    % Generate test signal:
    N   = 44100;
    fs  = 44100;
    lvl_dbfs = -3;
    sig1 = wgn(2*N,1,lvl_dbfs); 
    t = ( -length(sig1)/2:length(sig1)/2-1 )/fs;
    
    %% Interference 2 sources with a separation of d at mic 1 (only freqs):
    theta     = atan( distS1_S2/distS1_MicD ); 
    theta_deg = theta/pi*180; 
    
    d       = distS1_S2;
    n       = [1 2 3];
    deltat  = d*sin(theta)/c;
    
    M = round(deltat * fs);
    sig2 = [zeros(M-1,1); sig1(1:length(sig1)-M+1)];
    
    idx = find(t>=0);

    y1 = sig1(idx)+sig2(idx); % no alpha, since no floor is involved

    info.fs       = fs;
    info.typeplot = 3;
    info.bNewFigure = 1;
    freqfft(y1,K,info);
    
    hold on
    % fmax = n*c/(d*sin(theta));
    % % or
    fmax =    n   /   deltat ;
    fmin = (2*n-1)/(2*deltat);
    
    for i = 1:length(fmin)
        plot([fmin(i) fmin(i)], [15 50],'r--')
        plot([fmax(i) fmax(i)], [15 50],'ro-')
    end
    
    %% Interference 1 source with its first reflection at mic 1 (white noise):
    
    deltat = (distS2im_MicD - distS2_MicD)/c;

    M = round(deltat * fs);
    sig2 = [zeros(M-1,1); sig1(1:length(sig1)-M+1)];
    
    idx = find(t>=0);

    y2 = sig1(idx)+alfa*sig2(idx);

    figure;
    subplot(3,1,1)
    info.fs       = fs;
    info.typeplot = 3;
    info.bNewFigure = 0;
    freqfft(y2,K,info);
    xlim([20 fmaxPlot])
    hold on;
    
    title('S_2 on distant mic')
    
    n2 = [1:10];
    fmax2 = n2/deltat;
    fmin2 = (2*n2-1)/(2*deltat);
    
    deltaf2 = fmin2(2)-fmin2(1);
    
    for i = 1:length(fmin2)
        plot([fmin2(i) fmin2(i)], [15 50],'r--')
        plot([fmax2(i) fmax2(i)], [15 50],'ro-')
    end
    
    %% Interference 1 source with its first reflection at mic 1 (white noise):
    
    deltat = (distS1im_MicD - distS1_MicD)/c;

    M = round(deltat * fs);
    sig3 = [zeros(M-1,1); sig1(1:length(sig1)-M+1)];
    
    idx = find(t>=0);

    y3 = sig1(idx)+alfa*sig3(idx);

    subplot(3,1,2)
    info.fs       = fs;
    info.typeplot = 3;
    info.bNewFigure = 0;
    freqfft(y3,K,info);
    xlim([20 fmaxPlot])
    hold on;
    
    title('S_1 on distant mic')
    
    n2 = [1:10];
    fmax3 = n2/deltat;
    fmin3 = (2*n2-1)/(2*deltat);
    deltaf3 = fmin3(2)-fmin3(1);
    
    for i = 1:length(fmin2)
        plot([fmin3(i) fmin3(i)], [15 50],'r--')
        plot([fmax3(i) fmax3(i)], [15 50],'ro-')
    end
    
    %% Sum of all interfered signals:
    
    yt = (-2*sig1(idx) + y1 + y2 + y3)/3;
    
    subplot(3,1,3)
    info.fs       = fs;
    info.typeplot = 3;
    info.bNewFigure = 0;
    freqfft(yt,K,info);
    xlim([20 fmaxPlot])
    hold on;
    
    disp('')
    
end

% Not tested Unix:
bDoCalibration = 0;

%% Calibration
if bDoCalibration
    
    % dir_audio = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\03-Wav-files-1-referentie\';
    dir_audio = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\1-referentie-BP-filtered\at441kHz\';
    dir_audio_mod_ane   = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data\at44100Hz\';
    dir_audio_modelled  = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-reflections\at44100Hz\';
    output_dir  = [Get_TUe_paths('outputs') 'Hummer-audio-files' delim];
    
    field       = [1 2]; % 1 = distant mic; 2 = close mic 
    field_lbl   = {'dist','close'};
    
    ac_mode     = [2 3 5];
    ac_mode_idx = 1:length(ac_mode);
    
    cal_level   = [60 75 82]; % ac mode 2, 3, 5
    delta_dB_far = [-6 0]; % distance distant mic = 2 * close mic
    
    ver_anechoic = [2 1 1]; % ac-mode-5, take 1 = contaminated with ac mode 4
    ver_reverberant = [1 2 3]; % ac-mode-5
    
    freqs2cal   = [ 400 400; ... 
                    630 630; ...
                    1000 1250];
    
    %% Measurements:
    % Anechoic conditions:
    for i = field
        
        for j = ac_mode_idx
            
            files = sprintf('modus-%.0f_v%.0f-%.0f-filt-new.wav',ac_mode(j)-1,ver_anechoic(ac_mode_idx(j)), i); % e.g. modus-4_v3-1filt-new.wav = ac mode 5, distant mic
            fullfile = [dir_audio files];
            
            tmp.nAnalyser = 10; % One-third-OB
            tmp.CalMethod = 1; % AMT, 100 dB = 0 dBFS
            tmp.bPlot = 0;
            out = PsySoundCL(fullfile,tmp);
            f   = out.f;
            idx = find( f >= freqs2cal(j,1) & f <= freqs2cal(j,2));
            current_dB = sum_db( out.DataSpecOneThirdAvg(idx) );
            corr_factor_dB(j) = cal_level(j) - current_dB + delta_dB_far(i);
            corr_factor(j)     = From_dB(corr_factor_dB(j));
            
            [x fs] = Wavread([dir_audio files]);
            y2 = corr_factor(j)*x;
            
            fileout = sprintf('%smeas-ac-mode-%.0f-%s-anechoic', output_dir,ac_mode(j),field_lbl{i});
            
            Wavwrite(y2, fs, fileout)
            
            disp('')
            
        end
        
    end
    
    % Reverberant conditions:
    for i = field
        
        for j = ac_mode_idx
            
            files = sprintf('modus-%.0f_v%.0f-%.0f-filt-new.wav',ac_mode(j)-1,ver_reverberant(ac_mode_idx(j)), i); % e.g. modus-4_v3-1filt-new.wav = ac mode 5, distant mic
            fullfile = [dir_audio files];
            
            tmp.nAnalyser = 10; % One-third-OB
            tmp.CalMethod = 1; % AMT, 100 dB = 0 dBFS
            tmp.bPlot = 0;
            out = PsySoundCL(fullfile,tmp);
            f   = out.f;
            idx = find( f >= freqs2cal(j,1) & f <= freqs2cal(j,2));
            current_dB = sum_db( out.DataSpecOneThirdAvg(idx) );
            corr_factor_dB(j) = cal_level(j) - current_dB;
            corr_factor(j)     = From_dB(corr_factor_dB(j));
            
            [x fs] = Wavread([dir_audio files]);
            y2 = corr_factor(j)*x;
            
            fileout = sprintf('%smeas-ac-mode-%.0f-%s-reverberant', output_dir,ac_mode(j),field_lbl{i});
            
            Wavwrite(y2, fs, fileout)
            
            disp('')
            
        end
        
    end
    
    %% Predicted:
    for i = field
        
        for j = ac_mode_idx
            
            files = sprintf('modus-%.0f-v_%.0f-new.wav',ac_mode(j)-1, i); % e.g. modus-4_v3-1filt-new.wav = ac mode 5, distant mic
            fullfile = [dir_audio_mod_ane files];
            
            tmp.nAnalyser = 10; % One-third-OB
            tmp.CalMethod = 1; % AMT, 100 dB = 0 dBFS
            out = PsySoundCL(fullfile,tmp);
            f   = out.f;
            idx = find( f >= freqs2cal(1) & f <= freqs2cal(2));
            current_dB = sum_db( out.DataSpecOneThirdAvg(idx) );
            corr_factor_dB(j) = cal_level(j) - current_dB + delta_dB_far(i);
            corr_factor(j)     = From_dB(corr_factor_dB(j));
            
            [x fs] = Wavread([dir_audio_mod_ane files]);
            y2 = corr_factor(j)*x;
            
            fileout = sprintf('%smodel-ac-mode-%.0f-%s-anechoic', output_dir,ac_mode(j),field_lbl{i});
            
            Wavwrite(y2, fs, fileout)
            
            disp('')
            
        end
        
    end
    
    % Reverberant conditions:
    for i = field
        
        for j = ac_mode_idx
            
            files = sprintf('modus-%.0f-v_%.0f-new.wav',ac_mode(j)-1, i); 
            fullfile = [dir_audio_modelled files];
            
            tmp.nAnalyser = 10; % One-third-OB
            tmp.CalMethod = 1; % AMT, 100 dB = 0 dBFS
            out = PsySoundCL(fullfile,tmp);
            f   = out.f;
            idx = find( f >= freqs2cal(1) & f <= freqs2cal(2));
            current_dB = sum_db( out.DataSpecOneThirdAvg(idx) );
            corr_factor_dB(j) = cal_level(j) - current_dB;
            corr_factor(j)     = From_dB(corr_factor_dB(j));
            
            [x fs] = Wavread([dir_audio_modelled files]);
            y2 = corr_factor(j)*x;
            
            fileout = sprintf('%smodel-ac-mode-%.0f-%s-reverberant', output_dir,ac_mode(j),field_lbl{i});
            
            Wavwrite(y2, fs, fileout)
            
            disp('')
            
        end
        
    end
    
end


if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
