function [outsig Fi Ai sigma_i Phi_i L Ai_dB] = Get_ESPRIT_piano(insig,fs,abovedB,M,L)
% function [outsig Fi Ai sigma_i Phi_i L Ai_dB] = Get_ESPRIT_piano(insig,fs,abovedB,M,L)
%
% 1. Description:
%       
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 07/05/2016
% Last update on: 07/05/2016 
% Last use on   : 07/05/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    abovedB = -100;
end

if nargin < 4
    M = min(2^12,length(insig));
end

if nargin < 5
    L = 200;
end

p = 2*L;
[outsig Fi Ai sigma_i Phi_i L] = Get_ESPRIT_analysis(insig,p,M,fs);

idxmax = find(Fi<10000,1,'last'); % fmax is never a value above 10 kHz
Ai_dB = To_dB(Ai) - max(To_dB(Ai(1:idxmax)));

idx = find(Ai_dB < abovedB);
Fi(idx) = [];
Ai(idx) = [];
sigma_i(idx) = [];
Phi_i(idx) = [];
Ai_dB(idx) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
