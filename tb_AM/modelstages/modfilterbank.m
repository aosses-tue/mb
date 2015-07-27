function [outsig,mfc] = modfilterbank(insig,fs,fc,varargin)
% function [outsig,mfc] = modfilterbank(insig,fs,fc,varargin)
% 
% 1. Descrption: 
%   Modulation filter bank
% 
%   Usage: [outsig, mfc] = modfilterbank(insig,fs,fc);
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
%   MODFILTERBANK(insig,fs,fc) applies a modulation filterbank to the input
%   signals insig which are sampled with a frequency of fs Hz. Each column in
%   insig is assumed to be bandpass filtered with a center frequency stored in fc.
%
%   By default, the modulation filters will have center frequencies
%   0,5,10,16.6,27.77,... where each next center frequency is 5/3 times the
%   previous one. For modulation frequencies below (and including) 10 Hz,
%   the real value of the filters are returned, and for higher modulation
%   center frequencies, the absolute value (the envelope) is returned.
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

% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
% AUTHOR: Stephan Ewert
%
% Modifications by Morten L. Jepsen and Peter L. Søndergaard.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

definput.keyvals.mfc=[];
[flags,kv]=ltfatarghelper({},definput,varargin);

nfreqchannels   = length(fc); % number of frequency channels of the auditory filterbank
startmf         = 5; % fcm min = 5 Hz

Q = 2;
bw = 5;
ex=(1+1/(2*Q))/(1-1/(2*Q));

outsig=cell(nfreqchannels,1);

% second order modulation Butterworth LPF with a cut-off frequency of 2.5 Hz.
[b_lowpass,a_lowpass] = butter(2,2.5/(fs/2));

% first order modulation Butterworth LPF with a cut-off frequency of 150 Hz. 
% This is to remove all modulation frequencies above 150 Hz.
[b_highest,a_highest] = butter(1,150/(fs/2));

% Set the highest modulation frequency as proportion of the corresponding
% center frequency.
umf = min(fc.*0.25, 1000);  

for freqchannel=1:nfreqchannels

    % Cut away highest modulation frequencies
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
            % figure; t_tmp = (1:length(outtmp))/fs; plot(t_tmp,outtmp,t_tmp,outsigblock(:,nmfc))
            % legend('sig_i','sig_i filt'); grid on; title( sprintf('fc=%.1f [Hz]',mfc(nmfc)) )
        end
    
    end
  
  %% ------------ post-processing --------------------
  
  for nmfc=1:length(mfc) % v2 MJ 17. oct 2006
    if mfc(nmfc) <= 10 % f below 10 Hz
      outsigblock(:,nmfc) = 1*real(outsigblock(:,nmfc));
    else
      outsigblock(:,nmfc) = 1/sqrt(2)*abs(outsigblock(:,nmfc));
    end
  end
  
  outsig{freqchannel}=outsigblock;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ------------ subfunctions ------------------------

% complex frequency shifted first order lowpass
function [b,a] = efilt(w0,bw);

e0 = exp(-bw/2);

b = 1 - e0;
a = [1, -e0*exp(1i*w0)];
