function r20141024_determine_thres(nStimuli)
% function r20141024_determine_thres(nStimuli)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       nStimuli = 1; % determine thresholds for 10 ms test stimuli
%       nStimuli = 2; % determine thresholds for  5 ms test stimuli
%       r20141024_running_noise(nStimuli);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/10/2014
% Last update on: 27/10/2014 % Update this date manually
% Last use on   : 27/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

% method = 'dau1996'; % with overshoot limiting
method = 'dau1996'; % without overshoot limiting

%% Reading and adjusting onset of test stimuli

lvl_ref = 85; 

pathaudio = Get_TUe_paths('outputs');
% pathaudio = 'D:\Output\';

if nargin == 0
    close all
    nStimuli = 1;
end

switch nStimuli
    case 1
        % test signals with onset at 100 ms
        % noise        with onset at   0 ms 
        lvl_dB  = [66 70 78 82];
        filename_N      = [ pathaudio 'dau1996b_expI1_noisemasker.wav'];
        filename_S1   = [ pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_dB(1)) '.wav'];
        filename_S2   = [ pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_dB(2)) '.wav'];
        filename_S3   = [ pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_dB(3)) '.wav'];
        filename_S4   = [ pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_dB(4)) '.wav'];
        filename_Sref = [ pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_ref) '.wav'];
        onsetS        = [125]*1e-3; % fixed
        durS          = 5e-3;
        
    case 10 % check this experiment number
        lvl_dB  = [56 66 76 85];
        filename_N      = [ pathaudio 'dau1996b_expIB0_noisemasker.wav'];
        filename_S1   = [ pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_dB(1)) '.wav'];
        filename_S2   = [ pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_dB(2)) '.wav'];
        filename_S3   = [ pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_dB(3)) '.wav'];
        filename_S4   = [ pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_dB(4)) '.wav'];
        filename_Sref = [ pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_ref) '.wav'];
        onsetS        = [10]*1e-3;
        durS          = 10e-3;
    
end

[insig_N  fs] = Wavread(filename_N);
[insig_S1 fs] = Wavread(filename_S1);
[insig_S2 fs] = Wavread(filename_S2);
[insig_S3 fs] = Wavread(filename_S3);
[insig_S4 fs] = Wavread(filename_S4);

[insig_Sref fs] = Wavread(filename_Sref);

switch nStimuli
    case 1
        toPad = Gen_silence(0,fs); % No padding needed
        Ni    = round(onsetS * fs);
        
    case 10 % Check
        toPad = Gen_silence(onsetS,fs);
        Ni    = 1;
    
end

N_added = length(toPad);

insig1  = [toPad; insig_S1(1:end-N_added)];
insig2  = [toPad; insig_S2(1:end-N_added)];
insig3  = [toPad; insig_S3(1:end-N_added)];
insig4  = [toPad; insig_S4(1:end-N_added)];
insigref  = [toPad; insig_Sref(1:end-N_added)];
        
durN = 200e-3;
hopsize = 10e-3*fs;
bufsize = 20e-3*fs;

%%

[RM  fc t] = Get_internal_representations(insig_N          ,fs,method);
[RMT     ] = Get_internal_representations(insig_N+insigref ,fs,method);

[RMTc1   ] = Get_internal_representations(insig_N+insig1 ,fs,method);
[RMTc2   ] = Get_internal_representations(insig_N+insig2 ,fs,method);
[RMTc3   ] = Get_internal_representations(insig_N+insig3 ,fs,method);
[RMTc4   ] = Get_internal_representations(insig_N+insig4 ,fs,method);


idx = 13;
template = Get_template(RM,RMT,'template');

d(1) = Get_template(RM,RMTc1,'dprime',idx);
d(2) = Get_template(RM,RMTc2,'dprime',idx);
d(3) = Get_template(RM,RMTc3,'dprime',idx);
d(4) = Get_template(RM,RMTc4,'dprime',idx);
d = abs(d);

cc(1) = Get_template(template,RMTc1,'cross-correlation',idx);
cc(2) = Get_template(template,RMTc2,'cross-correlation',idx);
cc(3) = Get_template(template,RMTc3,'cross-correlation',idx);
cc(4) = Get_template(template,RMTc4,'cross-correlation',idx);

critD   = 1.26;
thrD    = interp1(d, lvl_dB, critD);

critCC = 6.5;
thrCC   = interp1(cc, lvl_dB, critCC);

%% Figures

% Noise, waveform
N = round(durN*fs);
figure;
plot(t(1:N),insig_N(1:N))
title(  sprintf( 'rms=%.2f dB',rmsdb(insig_N(Ni:Ni+N))+100 )  );

N = round(durS*fs);

% Test signals, wavefom
figure;
subplot(4,1,1)
plot(t(Ni:Ni+N),insig_S1(Ni:Ni+N))
title(  sprintf( 'rms=%.2f dB',rmsdb(insig_S1(Ni:Ni+N))+100 )  );

subplot(4,1,2)
plot(t(Ni:Ni+N),insig_S2(Ni:Ni+N))
title(  sprintf( 'rms=%.2f dB',rmsdb(insig_S2(Ni:Ni+N))+100 )  );

subplot(4,1,3)
plot(t(Ni:Ni+N),insig_S3(Ni:Ni+N))
title(  sprintf( 'rms=%.2f dB',rmsdb(insig_S3(Ni:Ni+N))+100 )  );

subplot(4,1,4)
plot(t(Ni:Ni+N),insig_S4(Ni:Ni+N))
title(  sprintf( 'rms=%.2f dB',rmsdb(insig_S4(Ni:Ni+N))+100 )  );

% Template
figure;
plot(t,template(:,idx))
title(sprintf('template, fc=%.1f',fc(idx)))

% Internal representations:
figure;
subplot(4,1,1)
plot(t,RM(:,idx),t,RMTc1(:,idx))
title(sprintf('IR, fc=%.1f,lvl=%.1f',fc(idx),lvl_dB(1)))

subplot(4,1,2)
plot(t,RM(:,idx),t,RMTc2(:,idx))
title(sprintf('IR, fc=%.1f,lvl=%.1f',fc(idx),lvl_dB(2)))

subplot(4,1,3)
plot(t,RM(:,idx),t,RMTc3(:,idx))
title(sprintf('IR, fc=%.1f,lvl=%.1f',fc(idx),lvl_dB(3)))

subplot(4,1,4)
plot(t,RM(:,idx),t,RMTc4(:,idx))
title(sprintf('IR, fc=%.1f,lvl=%.1f',fc(idx),lvl_dB(4)))

% [statsN tbuf]   = Get_stats_from_audio_excerpt(RM(:,idx) , fs, bufsize, hopsize);

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
