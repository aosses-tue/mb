function [outsig,mfc] = modfilterbank1997(insig,fs,fc,varargin)
% function [outsig,mfc] = modfilterbank1997(insig,fs,fc,varargin)
% 
% 1. Descrption: 
%       Modulation filter bank as used in Dau1997.
% 
%   Input parameters:
%      insig  : Input signal(s)
%      fs     : Sampling rate in Hz,
%      fc     : Center frequencies of the input signals
%
%   Output parameters:
%      outsig : Modulation filtered signals
%      mfc    : Center frequencies of the modulation filters.
%
%   MODFILTERBANK1997(insig,fs,fc) applies a modulation filterbank to the 
%       input signals insig which are sampled with a frequency of fs Hz. 
%       Each column in insig is assumed to be bandpass filtered with a 
%       center frequency stored in fc.
%
%       By default, the modulation filters will have center frequencies 0 
%       (cut-off at 2.5 Hz),5,10,16.6,27.77,... where each next center 
%       frequency is 5/3 times the previous one. For modulation frequencies 
%       below (and including) 10 Hz, the real value of the filters are 
%       returned.
%  
%   References:
%     T. Dau, B. Kollmeier, and A. Kohlrausch. Modeling auditory processing
%     of amplitude modulation. I. Detection and masking with narrow-band
%     carriers. J. Acoust. Soc. Am., 102:2892-2905, 1997a.
%     
%     R. Fassel and D. Pueschel. Modulation detection and masking using
%     deterministic and random maskers. Contributions to Psychological
%     Acoustics, edited by A. Schick (Universitaetsgesellschaft Oldenburg,
%     Oldenburg), pages 419-429, 1993.
%
%   See also: dau1997preproc
%
%   Url: http://amtoolbox.sourceforge.net/doc/modelstages/modfilterbank.php
%
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
%  
% AUTHOR: Stephan Ewert
%
% Modifications by Morten L. Jepsen and Peter L. Soendergaard.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.keyvals.mfc=[];
[flags,kv]=ltfatarghelper({},definput,varargin);

nfreqchannels   = length(fc); % number of frequency channels of the auditory filterbank
startmf         = 5; % fcm min = 5 Hz

Q   = 2;
bw  = 5;
ex  = (1+1/(2*Q))/(1-1/(2*Q));

outsig=cell(nfreqchannels,1);

% second order modulation Butterworth LPF with a cut-off frequency of 2.5 Hz.
[b_lowpass,a_lowpass] = butter(2,2.5/(fs/2));
[b_highest,a_highest] = butter(1,150/(fs/2));
% % Set the highest modulation frequency as proportion of the corresponding
% % center frequency.
umf = 1000;  

for freqchannel=1:nfreqchannels

    %outtmp = insig(:,freqchannel);
    outtmp = filter(b_highest,a_highest,insig(:,freqchannel));
    if umf(freqchannel)==0
        % ----------- only lowpass ---------------------
        outsigblock = filter(b_lowpass,a_lowpass,outtmp);
        mfc = 0;

    else
        tmp = fix((min(umf(freqchannel),10) - startmf)/bw);
        tmp = 0:tmp;
        mfc = startmf + 5*tmp;
        tmp2 = (mfc(end)+bw/2)/(1-1/(2*Q));
        tmp = fix(log(umf(freqchannel)/tmp2)/log(ex));
        tmp = 0:tmp;
        tmp = ex.^tmp;
        mfc=[0 mfc tmp2*tmp];

        % --------- lowpass and modulation filter(s) ---
        outsigblock = zeros(length(insig),length(mfc)); % memory allocation
        outsigblock(:,1) = filter(b_lowpass,a_lowpass,outtmp);

        for nmfc=2:length(mfc)
            w0 = 2*pi*mfc(nmfc)/fs;
            if mfc(nmfc) < 10   % frequencies below 10 Hz.  - 8 Hz LPF to preserve modulation phase
                [b3,a3] = efilt(w0,2*pi*bw/fs);
            else
                % frequencies above 10 Hz: 
                %   - Dau1997b, pp.2894 [...] Within the model only the (Hilbert) envelope
                %     of the modulation filter outputs for center frequencies above
                %     10 Hz is further examined. [...]
                [b3,a3] = efilt(w0,w0/Q);
            end
      
            outsigblock(:,nmfc) = 2*filter(b3,a3,outtmp);
            
        end
    
    end
  
  %% ------------ post-processing --------------------
  
    for nmfc=1:length(mfc)
        if mfc(nmfc) <= 10 % f below 10 Hz
            outsigblock(:,nmfc) = 1*real(outsigblock(:,nmfc));
        else
            outsigblock(:,nmfc) = abs(outsigblock(:,nmfc)); 
        end
    end
    outsig{freqchannel}=outsigblock;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inline functions:

% 1. Complex frequency shifted first order low-pass
function [b,a] = efilt(w0,bw);

e0 = exp(-bw/2);

b = 1 - e0;
a = [1, -e0*exp(1i*w0)];
