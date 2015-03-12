function [modfcs outMTF SNRi m] = MTF(cleanEnvelope, processedEnvelope,fs)
% function [modfcs outMTF SNRi m] = MTF(cleanEnvelope, processedEnvelope,fs)
% 
% This function calculates the MTF in third octave bands of a signal
% envelope as described in Goldsworthy and Greenberg (2004)
%
%  processedEnvelope : The envelope of a signal which have been processed in
%                   some way
%  cleanEnvelope     : The envelope of the unproccesed signal
%  outMFT            :  vector with the outputs MTF as function of modulation
%  centerfrequency
%  Hx_oct            :  matrix with the envelope spectra of the input
%  stimuli
%  
%  SNRi       :  vector of aSNRs in dB, one for each center frequency as
%  Eqn. (2) in goldsworthy
%  Søren Jørgensen, October 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Center frequencies of the modulation filters
modfcs = [0.63 .8 1 1.25 1.6 2.0 2.5 3.15 4.0 5.0 6.3 8 10 12.5];

N = length(cleanEnvelope);

nstim = 2;
stim = zeros(N,nstim);
stim(:,1)= cleanEnvelope;
stim(:,2) = processedEnvelope;

% pads = zeros(1,N*3);
% N2 = N+length(pads);
% Ly2=pow2(nextpow2(N2));
% N = N2;

% Performing a third-octave band analysis:
for k = 1:nstim 
    [tmp(:,k)] =  ThirdOctave_rmsAnalysis(stim(:,k),fs,modfcs); %rms-levels in 1/3 octave filters
    % normalizing with mean intensity of  signal envelope to get modulation index values
    m(:,k) = tmp(:,k)/(mean(stim(:,k)));
end

% Calculation of MTF:
a = 1;
outMTF = a * ( m(:,2) ./m(:,1));
outMTF = min(outMTF,1);

% calculation of apparentSNR
SNRi = 10*log10(((outMTF) ./(1- outMTF)));
% truncation between +- 15 dB:
SNRi = min(SNRi,15);
SNRi = max(SNRi,-15);

% figure
% plot(stim(1,:)),hold on
% plot(stim(2,:),':r')
% % plot(stim(3,:),'.- g')
% plot(stim(4,:),'c')
% plot(noiseFloor,'black :')
% std(stim(3,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
