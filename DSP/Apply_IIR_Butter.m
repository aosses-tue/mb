function outsig = Apply_IIR_Butter(insig,fs,fc,type,order)
% function outsig = Apply_IIR_Butter(insig,fs,type,fc,order)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 18/03/2016
% Last update on: 18/03/2016 
% Last use on   : 18/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    type = 'low';
end

if nargin < 5
    order = 4; % slope order*6 dB/Oct
end

wc = fc/(fs/2); % normalised frequency (fs/2 = 1)
[b, a] = butter(order,wc,type);

outsig = filtfilt(b,a,insig);

if nargout == 0
    
    N = 8192;
    [H   ,F]=freqz(b,a,N,fs);
    [Hin ,F]=freqz(insig,1,N,fs);
    [Hout,F]=freqz(outsig,1,N,fs);
    
    phi     = phasez(b,a,N,fs);
    phi_in  = phasez( insig,1,N,fs);
    phi_out = phasez(outsig,1,N,fs);
    
    figure
    subplot(2,1,1)
    semilogx(F,20*log10(abs(H))); hold on
    semilogx(F,20*log10(abs(Hin)),'r');
    semilogx(F,20*log10(abs(Hout)),'g');
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [dB]')
    ha = gca;
    
    subplot(2,1,2)
    semilogx(F,phi); hold on
    semilogx(F,phi_in ,'r');
    semilogx(F,phi_out,'g');
    xlabel('Frequency [Hz]')
    ylabel('Phase [rad]')
    ha(end+1) = gca;
    linkaxes(ha,'x');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
