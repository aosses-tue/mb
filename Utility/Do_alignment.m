function [fta ta] = Do_alignment(t, ft, ti)
% function [fta ta] = Do_alignment(t, ft, ti)
%
% 1. Description:
%   t - time
%   ft - f(t)
%   ti - actual initial time
% 
%   ta - time aligned
%   fta - ft aligned
% 
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%       See also: r20151223_update_Antoine.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 31/07/2014
% Last update on: 23/12/2015 
% Last use on   : 23/12/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_t = find( t >= ti ,1,'first'); % the same as old: min(find( t >= ti ));
fta   = ft(idx_t:end);
L     = length(fta);
 
fta   = Do_truncate(fta,L);
ta    = Do_truncate(t ,L);

if nargout == 0
    figure;
    
    plot(t,ft,'b', ta,fta,'r--');
    grid on
    legend('no align', 'align')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
