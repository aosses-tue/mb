function [Ns Ns_MA] = Reduce_Bark_resolution(res)
% function [Ns Ns_MA] = Reduce_Bark_resolution(res)
%
% 1. Description:
%       Applies moving average filter to Specific loudness (Zwicker's model)
%       and then reduces frequency resolution (default from 0.1 Barks to 1
%       1 Bark)
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 3/7/2014
% Last update on: 3/7/2014 % Update this date manually
% Last used on  : 3/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns_MA = filter(res.b, res.a, res.InstantaneousSpecificLoudness); % Applying moving average

Step = 1/res.BarkStep;
tmpf = [];

for i = 1 : Step : size(Ns_MA,2)
    
    tmp     = Ns_MA(:,i:i+Step-1);
    tmpf    = [tmpf; sum(tmp')*res.BarkStep];
    
end

Ns = tmpf';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end