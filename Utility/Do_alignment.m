function [fta ta] = Do_alignment(t, ft, ti)
% function [y taligned] = Do_alignment(t, ft, ti)
%
% 1. Description:
%   t - time
%   ft - f(t)
%   ti - actual initial time
% 
%   ta - time aligned
%   fta - ft aligned
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 31/7/2014
% Last update on: 31/7/2014 % Update this date manually
% Last use on   : 31/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_t  = min(find( t >= ti ));

fta = ft(idx_t:end);

L = length(fta);
 
fta     = Do_truncate(fta,L);
ta      = Do_truncate(t ,L);

if nargout == 0
    figure;
    
    plot(t,ft,'b', ta,fta,'r--');
    grid on
    legend('no align', 'align')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
