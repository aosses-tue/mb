function [b a] = Create_Butter_att_oct(fc, fs, Att_per_oct,type)
% function [b a] = Create_Butter_att_oct(fc, fs, Att_per_oct,type)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 14/12/2015
% Last update on: 14/12/2015 
% Last use on   : 14/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order=ceil(Att_per_oct/6); % butterworth: A = 6 dB/Oct * Order

if nargin < 4
    if length(fc) == 1
        type = 'low'; % or 'high'
    elseif length(fc) == 2
        type = 'band-pass';
    end 
end

if length(fc) == 1
    [b,a]=butter(order, fc/(fs/2),type);
    % [b,a]=butter(order, fc/fs*2,'high');
elseif length(fc) == 2
    [b,a]=butter(order/2, fc/fs*2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
