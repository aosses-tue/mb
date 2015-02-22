function r20150130_update
% function r20150130_update
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 27/01/2015
% Last update on: 02/02/2015 % Update this date manually
% Last use on   : 02/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 1;
Diary(mfilename,bDiary);

bGetT   = 0; % Part 1
bLoadFiles = 0; % Part 2
bLearningCombFilter = 0; % Part 3
bDoCalibration = 1;
% bSave   = 0; 

paths   = Get_TUe_subpaths('db_voice_of_dragon');

%% Part 1
if bGetT
   
    Get_VoD_params(bGetT,bSave);
    Get_VoD_params_ledenmaat(bSave);
    
end

%% Part 2
if bLoadFiles
    
    ha = [];
%     x1 = [];
%     dir_1_ref       = paths.dir_measurements_cal{1};
%     dir_6_maateraf  = paths.dir_measurements_cal{6};
   
    f2fm = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\03-Wav-files-1-referentie\modus-1_v2-1filt.wav';
    f2cm = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\03-Wav-files-1-referentie\modus-1_v2-2filt.wav'; 
    f5fm = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\03-Wav-files-1-referentie\modus-4_v3-1filt.wav';
    f5cm = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\03-Wav-files-1-referentie\modus-4_v3-2filt.wav';
           
    f2fp = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\03-Wav-files-1-referentie\modus-1-v_1filt.wav';
    f2cp = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\03-Wav-files-1-referentie\modus-1-v_2filt.wav'; 
    f5fp = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\03-Wav-files-1-referentie\modus-4-v_1filt.wav';
    f5cp = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\03-Wav-files-1-referentie\modus-4-v_2filt.wav'; 
    
    T = [0.6 0.4 0.3 0.27];
    
%     f1 = {  [dir_1_ref         'modus-1_v2-2filt-new.wav'], ...
%             [dir_1_ref         'modus-1_v2-1filt-new.wav'], ...
%             [dir_1_ref         'modus-3_v4-2filt-new.wav'], ...
%             [dir_1_ref         'modus-4_v3-2filt-new.wav'] };
%     f6 = [dir_6_maateraf    'meas-ac-mode-2-maat-eraf-close-2-filt.wav'];
     
    [x2fm fs] = Wavread(f2fm);
    [x2cm fs] = Wavread(f2cm);
    [x2cp fs] = Wavread(f2cp);
    [x2fp fs] = Wavread(f2fp);

    [x5fm fs] = Wavread(f5fm);
    [x5cm fs] = Wavread(f5cm);
    [x5cp fs] = Wavread(f5cp);
    [x5fp fs] = Wavread(f5fp);

    t = ( 1:length(x2fm) )/fs;
    
%     info.fs = fs;
%     K = (4096)/2; % 4096-point FFT 
%     freqfft(x5cp(1:K),K*2,info);
%     xlim([900 1200])
%     
%     freqfft([x5cp(1:K); zeros(K,1)],K*2,info);
%     xlim([900 1200])
    
    %% Zero padding
    
    idx = find(t<0.25*5);
    
    % plt.figure(1, figsize=(9.5, 6))
    M = length(idx); % 128;
    tX = t(idx);
    
    N1 = 128;
    N2 = N1*32;
    N3 = N1*32*8;
    % x = cos(2*pi*2/M*[1:M]) .* transpose( hanning(M) );
    x = x5cp(1:M); 
    
    figure;
    subplot(4,1,1)
    plot(tX, x, 'b') %, marker='x', lw=1.5)
    title(sprintf('x, M=%.0f',M))
    xlabel('Time [s]')
    % axis([-M/2,M/2-1,-1,1])
 
    mX = 20 * log10(abs(fft(x, N1)));
    subplot(4,1,2)
    plot([0:N1-1], mX); %, marker='x', color='r', lw=1.5)
    axis([0,N1/2-1,-20,max(mX)+1])
    title(sprintf('magnitude spectrum: mX1, N=%.0f',N1))

    mX = 20 * log10(abs(fft(x, N2)));
    subplot(4,1,3)
    plot([0:N2-1], mX);
    axis([0,N2/2-1,-20,max(mX)+1])
    title(sprintf('magnitude spectrum: mX2, N=%.0f',N2))

    mX = 20 * log10(abs(fft(x, N3)));
    subplot(4,1,4)
    plot([0:N3-1],mX) %,marker='x',color='r', lw=1.5)
    axis([0,N3/2-1,-20,max(mX)+1])
    title(sprintf('magnitude spectrum: mX3, N=%.0f',N3))

    % plt.tight_layout()
    % plt.savefig('zero-padding.png')
    % plt.show()

    
    [xc1 fsc] = Wavread('D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\1-referentie-BP-filtered\modus 1_v1-2-filt.wav'); % starts at 0.2535 s
    xc2 = Wavread('D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\1-referentie-BP-filtered\modus 1_v2-2-filt.wav');
    tc = (1:length(xc1))/fsc;
    
    figure;
    plot(tc-0.27,xc1+0.05,tc,xc2-0.05)
    xlim([0 20-0.2535])
    
    figure;
    subplot(4,2,2)
    plot( t/T(1), x2cp ,'r');
    ha(end+1) = gca;
    title('referentie, close, ac mode 2, model')
    
    subplot(4,2,4)
    plot( t/T(1), x2fp ,'r');
    ha(end+1) = gca;
    title('referentie, distant, ac mode 2, model')
    
    subplot(4,2,1)
    plot( t/T(1), x2cm );
    ha(end+1) = gca;
    title('referentie, close, ac mode 2')
    
    subplot(4,2,3)
    plot( t/T(1), x2fm );
    ha(end+1) = gca;
    title('referentie, distant, ac mode 2')
    
    %
    subplot(4,2,6)
    plot( t/T(4), x5cp ,'r');
    ha(end+1) = gca;
    title('referentie, close, ac mode 5, model')
    
    subplot(4,2,8)
    plot( t/T(4), x5fp ,'r');
    ha(end+1) = gca;
    title('referentie, distant, ac mode 5, model')
    
    subplot(4,2,5)
    plot( t/T(4), x5cm );
    ha(end+1) = gca;
    title('referentie, close, ac mode 5')
    
    subplot(4,2,7)
    plot( t/T(4), x5fm );
    ha(end+1) = gca;
    title('referentie, distant, ac mode 5')
    
    linkaxes(ha,'x');
     xlim([0 2])
end

% Part 3
if bLearningCombFilter
    
    S1 = [0 0 2.23];
    S2 = [0,.67,2.23];
    MicD = [1.58,0,1.68];
    
    distS1_S2 = norm(S2-S1);
    distS1_MicD = norm(S1-MicD);
    distS2_MicD = norm(S2-MicD);
    theta = atan( distS1_S2/distS1_MicD )/pi*180;
    
    d = distS1_S2;
    
    dir ='D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2015-01-30-update\Hummer-audio-files\' 
    fi = [dir 'meas-ac-mode-2-close-anechoic.wav'];

    c = 340;
    
    [x fs] = Wavread(fi);
    % x = zeros(size(x));
    % x(1) = 0.8;
    d = 4.2/4;
    % d = 16.8;
    deltat = d/c;

    t = ( 1:length(x) )/fs;

    idx = find(t>deltat);
    x2 = [zeros(idx(1)-1,1); x(1:length(idx))];
    alfa = 1;

    y = x(1:length(x2))-alfa*x2;

    figure;
    plot(x+0.05), hold on
    plot(y,'r'), 

    Wavwrite(y,fs,sprintf('%stest-%.0f-0',dir,alfa*100))
end

%% Calibration
if bDoCalibration
    
    % dir_audio = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\03-Wav-files-1-referentie\';
    dir_audio = 'D:\Databases\dir01-Instruments\Voice-of-dragon\02-Wav-files\1-referentie-BP-filtered\at441kHz\';
    dir_audio_mod_ane   = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data\at44100Hz\';
    dir_audio_modelled  = 'D:\Databases\dir01-Instruments\Voice-of-dragon\03-Wav-files-predicted\01-Model\Data-new-reflections\at44100Hz\';
    output_dir  = [Get_TUe_paths('outputs') 'Hummer-audio-files' delim];
    
    field       = [1 2]; % 1 = distant mic; 2 = close mic 
    field_lbl   = {'dist','close'};
    
    ac_mode     = [2 3 4 5];
    ac_mode_idx = 1:length(ac_mode);
    
    cal_level   = [60 75 78 82]; % ac mode 2, 3, 5
    delta_dB_far = [-6 0]; % distance distant mic = 2 * close mic
    
    ver_anechoic    = [2 1 3 1]; % ac-mode-5, take 1 = contaminated with ac mode 4
    ver_reverberant = [1 2 1 3]; % ac-mode-5
    
    freqs2cal   = [400 1250];
    
    %% Measurements:
    % Anechoic conditions:
    for i = field
        
        for j = ac_mode_idx
            
            files = sprintf('modus-%.0f_v%.0f-%.0f-filt-new.wav',ac_mode(j)-1,ver_anechoic(ac_mode_idx(j)), i); % e.g. modus-4_v3-1filt-new.wav = ac mode 5, distant mic
            fullfile = [dir_audio files];
            
            tmp.nAnalyser = 10; % One-third-OB
            tmp.CalMethod = 1; % AMT, 100 dB = 0 dBFS
            out = PsySoundCL(fullfile,tmp);
            f   = out.f;
            idx = find( f >= freqs2cal(1) & f <= freqs2cal(2));
            current_dB = sum_db( out.DataSpecOneThirdAvg(idx) );
            corr_factor_dB(j) = cal_level(j) - current_dB + delta_dB_far(i);
            corr_factor(j)     = From_dB(corr_factor_dB(j));
            
            [x fs] = Wavread([dir_audio files]);
            y = corr_factor(j)*x;
            
            fileout = sprintf('%smeas-ac-mode-%.0f-%s-anechoic', output_dir,ac_mode(j),field_lbl{i});
            
            Wavwrite(y, fs, fileout)
            
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
            out = PsySoundCL(fullfile,tmp);
            f   = out.f;
            idx = find( f >= freqs2cal(1) & f <= freqs2cal(2));
            current_dB = sum_db( out.DataSpecOneThirdAvg(idx) );
            corr_factor_dB(j) = cal_level(j) - current_dB;
            corr_factor(j)     = From_dB(corr_factor_dB(j));
            
            [x fs] = Wavread([dir_audio files]);
            y = corr_factor(j)*x;
            
            fileout = sprintf('%smeas-ac-mode-%.0f-%s-reverberant', output_dir,ac_mode(j),field_lbl{i});
            
            Wavwrite(y, fs, fileout)
            
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
            y = corr_factor(j)*x;
            
            fileout = sprintf('%smodel-ac-mode-%.0f-%s-anechoic', output_dir,ac_mode(j),field_lbl{i});
            
            Wavwrite(y, fs, fileout)
            
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
            y = corr_factor(j)*x;
            
            fileout = sprintf('%smodel-ac-mode-%.0f-%s-reverberant', output_dir,ac_mode(j),field_lbl{i});
            
            Wavwrite(y, fs, fileout)
            
            disp('')
            
        end
        
    end
    
end


if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
