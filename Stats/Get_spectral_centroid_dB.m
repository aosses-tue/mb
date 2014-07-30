function fc = Get_spectral_centroid_dB(f,ydB,dB_below)
% function fc = Get_spectral_centroid_dB(f,ydB,dB_below)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 30/7/2014
% Last update on: 30/7/2014 % Update this date manually
% Last used on  : 30/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    dB_below = 20;
end

Maxi_y = max(ydB);
idx = find(ydB > Maxi_y-dB_below);

delta_idx = diff(idx);

f = f(idx);
y = From_dB(ydB(idx));

idx_n_lobes = find(delta_idx ~= 1);
N = length(idx_n_lobes);
if N > 1
    
    ti      = 1;
    tf      = idx_n_lobes(1);
    n = 1;
    sum_i   = 0;
    for i = 1:N-1
        idx_tmp1 = ti;
        idx_tmp2 = tf;
        
        sum_tmp = sum( From_dB(ydB(idx_tmp1:idx_tmp2)) );
        if sum_tmp > sum_i
            sum_i = sum_tmp;
            idx_i = idx_tmp1:idx_tmp2;
        end
        
        ti = idx_n_lobes(n  );
        tf = idx_n_lobes(n+1);
        
    end
    idx = idx_i;
    warning('Non consecutive sections...')
    f = f(idx);
    y = From_dB(ydB(idx));
    
elseif N == 1
    if round(idx_n_lobes/length(idx)) ~= 0
        idx = 1:idx_n_lobes; 
    else
        idx = idx_n_lobes:length(idx);
    end
end

fc = Get_spectral_centroid(f,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
