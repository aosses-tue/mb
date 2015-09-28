function r20141024_determine_thres(nStimuli,pathaudio)
% function r20141024_determine_thres(nStimuli,pathaudio)
%
% 1. Description:
%
%       Audio files needed, nStimuli = 1:
%           - dau1996b_expI_noisemasker.wav
%           - dau1996b_expI1_stim-5ms-66.wav
%           - dau1996b_expI1_stim-5ms-70.wav
%           - dau1996b_expI1_stim-5ms-74.wav
%           - dau1996b_expI1_stim-5ms-78.wav
%           - dau1996b_expI1_stim-5ms-82.wav
% 
% 2. Stand-alone example:
%       nStimuli = 1; % determine thresholds for 10 ms test stimuli
%       nStimuli = 2; % determine thresholds for  5 ms test stimuli
%       r20141024_running_noise(nStimuli);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 23/10/2014
% Last update on: 27/10/2014 % Update this date manually
% Last use on   : 19/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

% method = 'dau1996'; % with overshoot limiting
method = 'dau1996a'; % without overshoot limiting

%% Reading and adjusting onset of test stimuli

bNew = 1;
bOld = ~bNew;

lvl_ref = 85; 

if bOld
    if nargin < 2
        pathaudio = Get_TUe_paths('outputs');
    end
end

if bNew
    
    p = Get_date;
    if nargin < 2
        if ~isunix
            pathaudio = [Get_TUe_paths('outputs') mfilename p.date4files delim];
            Mkdir(pathaudio);
        else
            pathaudio = [Get_TUe_paths('outputs') 'r20141024_update_dau_et_al' delim 'Stimuli' delim];
        end
    end
    
end

% pathaudio = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Masking\dau1996b\Stimuli\';
% warning('Temporally changed')

if nargin == 0
    close all
    nStimuli = 1;
end

switch nStimuli
    case 1
        % test signals with onset at 100 ms
        % noise        with onset at   0 ms 
        lvl_dB  = [66 70 78 82];
        onsetS  = [125]*1e-3; % fixed
        durS    = 5e-3;
        
        if bNew
            filename_N    = [pathaudio 'dau1996b_expI_noisemasker.wav'];
        end
        
        if bNew | bOld
            filename_S1   = [pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_dB(1)) '.wav'];
            filename_S2   = [pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_dB(2)) '.wav'];
            filename_S3   = [pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_dB(3)) '.wav'];
            filename_S4   = [pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_dB(4)) '.wav'];
            filename_Sref = [pathaudio 'dau1996b_expI1_stim-5ms-' num2str(lvl_ref) '.wav'];
        end
        
        if bOld
            filename_N    = [pathaudio 'dau1996b_expI1_noisemasker.wav'];
        end
        
    case 10 % check this experiment number
        lvl_dB  = [56 66 76 85];
        onsetS  = [10]*1e-3;
        durS    = 10e-3;
        
        if bNew
            filename_N    = [pathaudio 'dau1996b_expI_noisemasker.wav'];
        end
        
        if bOld
            filename_N    = [pathaudio 'dau1996b_expIB0_noisemasker.wav'];
            filename_S1   = [pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_dB(1)) '.wav'];
            filename_S2   = [pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_dB(2)) '.wav'];
            filename_S3   = [pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_dB(3)) '.wav'];
            filename_S4   = [pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_dB(4)) '.wav'];
            filename_Sref = [pathaudio 'dau1996b_expIB0_stim-10ms-' num2str(lvl_ref) '.wav'];
            
        end
    
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

insig_s = insig_N+insigref;

% opts.bAddNoise = 0; % As reported on 24/10/2014
opts.bAddNoise = 1; 
opts.sigma = 20;
[RM  fc t] = Get_internal_representations(insig_N ,fs,method,opts);
[RMT     ] = Get_internal_representations(insig_s ,fs,method,opts);

insig_N1 = insig_N+insig1;
insig_N2 = insig_N+insig2;
insig_N3 = insig_N+insig3;
insig_N4 = insig_N+insig4;

[RMTc1   ] = Get_internal_representations(insig_N1 ,fs,method,opts);
[RMTc2   ] = Get_internal_representations(insig_N2 ,fs,method,opts);
[RMTc3   ] = Get_internal_representations(insig_N3 ,fs,method,opts);
[RMTc4   ] = Get_internal_representations(insig_N4 ,fs,method,opts);

N = 50; % PDF-resolution
M = 100;
idx = 13;

[PDF1 yi1] = Probability_density_function(insig_N1,N);
[PDF2 yi2] = Probability_density_function(insig_N2,N);
[PDF3 yi3] = Probability_density_function(insig_N3,N);
[PDF4 yi4] = Probability_density_function(insig_N4,N);
[PDFs yis] = Probability_density_function(insig_s ,N);
[PDFn yin] = Probability_density_function(insig_N,N);

[PDFR1 yr1] = Probability_density_function(RMTc1(:,idx),M);
[PDFR2 yr2] = Probability_density_function(RMTc2(:,idx),M);
[PDFR3 yr3] = Probability_density_function(RMTc3(:,idx),M);
[PDFR4 yr4] = Probability_density_function(RMTc4(:,idx),M);
[PDFRs yrs] = Probability_density_function(RMT(:,idx)  ,M);
[PDFRn yrn] = Probability_density_function(RM(:,idx)   ,M);


[PDIFF1 yd1] = Probability_density_function(RMTc1(:,idx)-RM(:,idx),M);
[PDIFF2 yd2] = Probability_density_function(RMTc2(:,idx)-RM(:,idx),M);
[PDIFF3 yd3] = Probability_density_function(RMTc3(:,idx)-RM(:,idx),M);
[PDIFF4 yd4] = Probability_density_function(RMTc4(:,idx)-RM(:,idx),M);
[PDIFFs yds] = Probability_density_function(RMT(:,idx)-RM(:,idx)  ,M);
% [PDIFFn ydn] = Probability_density_function(RM(:,idx)   ,M);

template = Get_template(RM,RMT,'template');

figure;
subplot(2,1,1)
plot(yin,PDFn, yis,PDFs), hold on, grid on; legend('noise','sig1')
title('Waveforms')

% subplot(2,1,2)
% plot(yrn,PDFRn, yrs,PDFRs), hold on, grid on; %legend('noise','sig1')
% title('Internal representations')

subplot(2,1,2)
plot(yd1,PDIFF1, yds,PDIFFs), hold on, grid on; legend('sig1','sig-supra')
title('Internal representations')

opts = [];
opts.idx = idx;
opts.fs = fs;
d(1) = Get_decision_criterion(RM,RMTc1,'dprime',opts);
d(2) = Get_decision_criterion(RM,RMTc2,'dprime',opts);
d(3) = Get_decision_criterion(RM,RMTc3,'dprime',opts);
d(4) = Get_decision_criterion(RM,RMTc4,'dprime',opts);
d = abs(d);

cc(1) = Get_decision_criterion(template,RMTc1,'cross-correlation',opts);
cc(2) = Get_decision_criterion(template,RMTc2,'cross-correlation',opts);
cc(3) = Get_decision_criterion(template,RMTc3,'cross-correlation',opts);
cc(4) = Get_decision_criterion(template,RMTc4,'cross-correlation',opts);

critD   = 1.26;
thrD    = interp1(d, lvl_dB, critD);

critCC = 6.5; % As reported on 24/10/2014.
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
