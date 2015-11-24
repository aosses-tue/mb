function [env noise w] = Create_piano_noise(insig,fs)
% function [env noise w] = Create_piano_noise(insig,fs)
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

w = randn(11*44100,1); % 11-s noise
env = Get_envelope_piano(insig,fs);

% w = setdbspl(w,lvl);

noise = w(1:length(env)).*env;
% noise = w(1:length(env));

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
