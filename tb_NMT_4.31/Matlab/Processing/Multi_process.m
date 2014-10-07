function outs = Multi_process(maps, in)

% Multi_process: Process one input signal through multiple processes.
%
% outs = Multi_process(maps, in)
%
% maps:  A cell array of processing maps (parameter structs).
% in:    An input signal.
% outs:  A cell array of output signals.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:length(maps)
    outs{n} = Process(maps{n}, in);
end    
