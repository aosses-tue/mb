function Add_paths(misc)
% function Add_paths(misc)
%
% 1. Description:
%   Add to path all directories specified in 'misc'
% 
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%   misc = Get_TUe_paths;
%   Add_paths(misc);
%   
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on: 13/5/2014
% Last update: 13/5/2014 % Update this date manually
% Last used: 13/5/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = fieldnames(misc);

for i = 1:length(fields)
    try
        if strcmp( misc.(fields{i})(end) , delim)
            addpath(misc.(fields{i}))
            disp(['Added to path: ' fields{i}])
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])