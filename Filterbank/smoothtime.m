function ysmooth = smoothtime(y,N)
% function ysmooth = smoothtime(y,N)
%
% 1. Description:
%   It does apply an N-point moving average filter.
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 21/05/2014
% Last update: 21/05/2014 % Update this date manually
% Last used: 21/05/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 % info does not exist
    N = 11; % use an odd number
end

ysmooth = [];
K = length(y);

bTruncated = 0;

for i = 1:K
    
    Ndist   = (N-1)/2;
    window  = zeros(K,1);
    
    liminf = i-Ndist;
    limsup = i+Ndist;
    
    if liminf < 1
        liminf = 1;
        bTruncated = 1;
    end
    if limsup > length(y)
        limsup = length(y);
        bTruncated = 1;
    end
    
    window(liminf:limsup) = 1;

    for j = 1:size(y,2)
        ysmooth(i,j) = sum(y(:,j).*window)/N;
        if bTruncated == 1
            ysmooth(i,j) = 0;
        end
    end
    bTruncated = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF - ' mfilename])