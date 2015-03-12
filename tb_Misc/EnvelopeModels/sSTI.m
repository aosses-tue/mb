function res = sSTI(cleanSig,procSig,fs,plotBand)
% function res = sSTI(cleanSig,procSig,fs,plotBand)
% 
% 1. Description:
%       This function calculates the speech-based STI as described in 
%       Houtgast & Steeneken (1985)
% 
%   res = sSTI(cleanSig,procSig,fs,plotBand)
% 
%  Inputs:
%   cleanSig : clean speech signal
%   procSig  : processed or noisy speech signal
%   fs       : sampling frequency
%   plotBand : number between 0-7 indicating which audioband is considered
%               for plotting the MTF. 0 indicates no plotting. A value
%               between 1 and 7 will plot the modulation spectra of speech
%               and noisy speech and the MTF for the corresponding audio
%               band.
%
%  Output:
%   res.STI        : the STI value
%   res.MTF        : the MTF of the stimuli
%   res.aSNRi      : a matrix of the apparent SNR as a function of the audio and modulation band
%   res.m_values   : a matrix of the modulation-index values as s function of the audio and modulation band
%   res.modBands   : The centerfrequencies of the modulation bands used in the MTF computation
%   res.audBands   : The centerfrequencies of the audio bands used in STI computation
% 
% Copyright Soren Jorgensen created august 2010
% last update april 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    plotBand = 0;
end

AudioBands  =[125  250  500  1000 2000 4000 8000];
wk          = [0.13 0.14 0.11 0.12 0.19 0.17 0.14]; % weights

stim(:,1) = cleanSig;
stim(:,2) = procSig;

nstim   = 2;
Nbands  = 7;

%% Calculating filterbank time domain output everything is put through filterbank
for k = 1:nstim
    [output_time(:,:,k) output_spec(:,:,k)] = OctFilterBank(stim(:,k),fs,AudioBands);
end

% Calculating envelopes of temporal outputs
for p = 1:Nbands
            
    % lowpass filtering at 40 Hz
    [bb, aa] = butter(2, 40*2/fs);
    D = 1;
    fs2 = fs/D;
    
    for k = 1:2
            
        env(:,k) = abs(hilbert(output_time(:,p,k))); %  envelope
        tmp = filter(bb,aa,env(:,k));
        %down sampling
        tmp = resample(tmp,fs2,fs);
        env_int(:,k) = tmp.^2; %envelope intensity

    end
    
    %% MTF
    [modBands outMTF(:,p) aSNRi(:,p) m_values(:,:,p)] = MTF(env_int(:,1), env_int(:,2),fs2);
         
end
    
aSNRk           = mean(aSNRi,1);
aSNR            = sum(aSNRk.*wk);
sti             = (aSNR+15)/30;
res.STI         = sti;
res.MTF         = outMTF;
res.aSNRi       = aSNRi;
res.m_values    = m_values;
res.modBands    = modBands;
res.audBands    = AudioBands;

if isempty(find(m_values(:,1)==0,1)) ==0
    disp('There was a 0 on in Hx_oct')
end
    
if plotBand>0
    
    band = plotBand;
    fnts = 14;
    lw = 2;
    figure
    subplot(2,1,1)
    plot((m_values(:,:,band)))
    % ylim([0 1.5])
    xlim([0 16])
    title(['modulation spectra ' num2str(AudioBands(band)) ' Hz octave band'])
    legend('Clean speech','Noisy speech')
    set(gca,'xTick',1:length(modBands),'xTickLabel',modBands,'ytick',0:.2:1.7,'FontSize',fnts,'FontWeight','b');
    xlabel('Modulation Band fc (Hz)')
    ylabel('Modulation index')
    
    subplot(2,1,2)
    plot(outMTF(:,band))
    ylim([0 1])
    xlim([0 16])
    title(['MTF ' num2str(AudioBands(band)) ' Hz octave band' ])
    set(gca,'xTick',1:length(modBands),'xTickLabel',modBands,'ytick',0:.1:1,'FontSize',fnts,'FontWeight','b');
    xlabel('Modulation Band fc (Hz)')
    ylabel('Modulation reduction')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
