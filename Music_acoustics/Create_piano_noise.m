function [env noise w] = Create_piano_noise(insig,fs,method)
% function [env noise w] = Create_piano_noise(insig,fs,method)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 24/11/2015
% Last update on: 24/11/2015 
% Last use on   : 24/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    method = 1;
end

w = randn(11*44100,1); % 11-s noise
env = Get_envelope_piano(insig,fs);
Menv = max(env);

switch method
    case 1
        noise = w(1:length(env)).*env;
                
    case 2

        N = round(500e-3*fs); % 500 ms
        K = N/2;
        idx = find(env > 0.1*Menv);
        insigbuf = buffer(insig(idx),N,K,'nodelay');
        
        [Pin F] = Get_PSD_analysis_arrays(insigbuf,fs);
        
        B = fir2(2^11,F/22050,sqrt(Pin));       %% sqrt because Pin is energy
        y = filter(B,1,w);
        noise = y(1:length(insig)).*env;
end

lvl = rmsdb(insig)+100;
noise = setdbspl(noise,lvl);

if nargout == 0
    t = (1:length(insig))/fs;

    figure;
    plot(t,noise,'r'); hold on;grid on
    plot(t,insig);  hold on;grid on
    plot(t,env,'r');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
