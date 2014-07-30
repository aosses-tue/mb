function fc = Get_spectral_centroid(f,y)
% function fc = Get_spectral_centroid(f,y)
%
% 1. Description:
%       fc - spectral centroid
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%   x = (11:16);
%   % y = 3*ones(1,6);
%   % y(1) = 1;
%   y = 2*[1 3 3 2 1 0];
%   Get_spectral_centroid(x,y);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 30/07/2014
% Last update on: 30/07/2014 % Update this date manually
% Last use on   : 30/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(y,1) ~= size(f,1)
    y = tranpose(y);
end

ysum = cumsum(y)/sum(cumsum(y));
ysum = cumsum(ysum);
idx_min = max(find(ysum<0.5));
fc = interp1(ysum(idx_min:idx_min+1), f(idx_min:idx_min+1),0.5);

if nargout == 0
    figure;
    plot(f,y,'b-',[fc fc],[min(y) max(y)],'r-')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
