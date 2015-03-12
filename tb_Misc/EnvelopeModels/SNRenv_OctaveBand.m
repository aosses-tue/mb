function res = SNRenv_OctaveBand(NoisySpeech,NoiseAlone,fs,plotBand)
% function res = SNRenv_OctaveBand(NoisySpeech,NoiseAlone,fs,plotBand)
% 
% 1. Description:
%           This function computes the SNRenv in seven octave bands
% 
% 
%   SNRenv_OctaveBand(NoisySpeech,NoiseAlone,fs,plotBand)
% 
%  Inputs:
%   NoisySpeech : The (processed) noisy speech signal 
%   NoiseAlone  : The (processed) Noise signal
%   fs       : sampling frequency
%   plotBand :  A value between 1 and 7 will plot the modulation excitation pattern of noisy speech
%               and noise alone for the corresponding audio band. 0 indicates no plotting.
%
%  Output:
%   res.SNRenv     : the overall SNRenv value integrated across audiobands
%   res.sEPSM_ExPtns : [7x2x7] matrix with Modulation excitation patterns
%                       of the stimuli. The first dimension is the modulation band, the second
%                       specifies whether its noisy speech (1) or noise alone (2) and the third
%                       indicates audio bands.
%   res.fcs_sEPSM   : The centerfrequencies of the modulation bands used in
%                       sEPSM
%   res.audBands   : The centerfrequencies of the audio bands used in STI
%                       computation
% 
% Copyright Søren Jørgensen created august 2010
% last update april 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    plotBand = 0;
end

AudioBands=[ 125 250 500 1000 2000 4000 8000 ];

nstim = 2;

stim(:,1) = NoisySpeech;
stim(:,2) = NoiseAlone;

Nbands = 7;

% Calculating filterbank time domain output everything is put through filterbank
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
            env(:,k) = resample(tmp,fs2,fs);
       
        end
    
     %% ------------estimation of SNRenv in this audio channel
    [fcs_sEPSM outSNRenvs(:,p) sEPSM_ExPtns(:,:,p)] = SNRenv_v1(env,fs);
        
        
    end
    
   tmp = outSNRenvs;
    %integrating across modulation and audio frequency
   combined_mod_tmp = sqrt(sum(tmp.^2,1));
   combined_aud  = (sqrt(sum(combined_mod_tmp.^2)));
   res.SNRenv = 10*log10(combined_aud);
   res.sEPSM_ExPtns = sEPSM_ExPtns;
   res.fcs_sEPSM = fcs_sEPSM; 
   res.audBands = AudioBands;
   
   
if plotBand > 0
    
   figure
   plot(10*log10(sEPSM_ExPtns(:,:,plotBand)))
   ylim([-25 10])
   xlim([0 8])
   ylabel('Envelope power (dB)')
   title(['Modulation excitation pattern for the audio filter tuned to ' num2str(AudioBands(plotBand)) ' Hz'])
   xlabel('Modulation filter f_c')
   legend('Noisy speech','Noise alone')
   set(gca, 'xtick',1:7, 'xticklabel',fcs_sEPSM)
    
end


