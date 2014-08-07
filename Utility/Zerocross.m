function [z time_info] = Zerocross(a,option)
% function [z time_info] = Zerocross(a,option)
%
% 1. Description:
%       a - MIR object
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 07/08/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    option = [];
    option.per = 'Second';
    option.dir = 'One'; % one direction
end

if iscell(a)
    a = a{1};
end

d = get(a,'Data');
f = get(a,'Sampling');
t = get(a,'Time');
while iscell(t)
    t = t{1};
end

% time_info.tup
% time_info.tdown
v = cell(1,length(d));

for h = 1:length(d)
    v{h} = cell(1,length(d{h}));
    for i = 1:length(d{h})
        di = d{h}{i};
        nc = size(di,2);
        nf = size(di,3);
        nl = size(di,1);
        bChangeDir = di(2:end,:,:).*di(1:(end-1),:,:) < 0;
        bUpwards    = ( di(2:end,:,:).*di(1:(end-1),:,:) < 0 ) & ( di(2:end,:,:) - di(1:(end-1),:,:) >= 0);
        bDownwards  = ( di(2:end,:,:).*di(1:(end-1),:,:) < 0 ) & ( di(2:end,:,:) - di(1:(end-1),:,:) <  0);
        idx_up      = find(bUpwards == 1);
        idx_down    = find(bDownwards == 1);
        
        time_info.t     = t;
        time_info.tup   = t(idx_up);
        time_info.tup_idx = idx_up;
        time_info.tdown = t(idx_down);
        time_info.tdown_idx = idx_down;
        
        zc = sum( bChangeDir ) /nl;
        if strcmp(option.per,'Second')
            zc = zc*f{h};
        end
        if strcmp(option.dir,'One')
            zc = zc/2;
        end
        v{h}{i} = zc;
    end
end
z = mirscalar(a,'Data',v,'Title','Zero-crossing rate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
