function [y_dB Fc] = Filterbank_analysis(x,fs,type)
% function [y_dB Fc] = Filterbank_analysis(x,fs,type)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on: 2/6/2014
% Last update: 01/07/2014 % Update this date manually
% Last used: 01/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y       = [];
y_dB    = [];
params.f0max = 2000;
params.f0min = 80;

switch type
    case 0 % OB analysis
        BandsPerOctave = 1;
        [y_dB Fc] = OB_analysis(x,fs,BandsPerOctave,params);
        
    case 1 % One-third OB analysis
        BandsPerOctave = 3;
        [y_dB Fc] = OB_analysis(x,fs,BandsPerOctave,params);
        
    case 2 % Gammatone analysis
     
        for k = 1:size(x,2)
            [aux Fc]    = Gammatone_analysis(x(:,k),fs);
            y_dB        = [y_dB transpose(rmsdb(aux))];
        end
        
end

Fc = Fc(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end