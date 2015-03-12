
%%
% This function calculates the SNRenv for the envelope signals
%
%  Env             : The envelopes used for calculating SNRenv
%  fs              : Sampling frequency
%  outSNRenvs      : 1x7 vector of overall SNRenv's , one for each modulation filter center frequency
%  fcs_EPSM        : center-frequencies of the modulation filterbank
%  sEPSM_ExPtns    : Envelope power Excitation patterns for of the input
% 
% Last edit: feb 2013 Søren Jørgensen 
%%
function [fcs_EPSM outSNRenvs sEPSM_ExPtns] = SNRenv_v1(Env,fs)

nstim =size(Env,2);


% Calculation of envelope power in 7 modulation filters
ExcPtn = zeros(7,nstim);
for k = 1:nstim
    [fcs_EPSM ExcPtn(:,k)] =  modFbank(Env(:,k),0,fs);
end

Noise_ExcPtn = ExcPtn(:,nstim);
Mix_ExcPtn =  ExcPtn(:,nstim-1);


% NaN values are set to zero
Noise_ExcPtn(isnan(Noise_ExcPtn)) = 0;

% Noisefloor cannot exceed the mix, since they are assumed to exist at the same time 
Noise_ExcPtn = min(Mix_ExcPtn,Noise_ExcPtn); 

% The noisefloor restricted to minimum 0.01 reflecting and internal noise
% threshold
Noise_ExcPtn = max(Noise_ExcPtn,0.01);
Mix_ExcPtn = max(Mix_ExcPtn,0.01);
% calculation of SNRenv
outSNRenvs = (((Mix_ExcPtn-Noise_ExcPtn ) ./Noise_ExcPtn));

% SNRenv - values are truncated to minimum -30 dB.
outSNRenvs = max(0.001,outSNRenvs);

% Excitation patterns
if nstim ==3
    sEPSM_ExPtns = [ExcPtn(:,1) Mix_ExcPtn Noise_ExcPtn];
else
    sEPSM_ExPtns = [Mix_ExcPtn Noise_ExcPtn];
end

