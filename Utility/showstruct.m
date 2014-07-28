function showstruct(thestruct, level)
%function showstruct(thestruct, level)
%
% 1. Description:
%       Show the fields belonging to a struct
%
% 2. Additional info:
%       Tested cross-platform: No
%       level  = 0; test OK
%       level ~= 0; not tested yet (01/07/2014)
%
% 3. Stand-alone example:
%       thestruct.bSave = 1;
%       thestruct.filename = 'audio4tests.wav';
%       showstruct(thestruct);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Tom Francart, ExpORL, KU Leuven, Belgium
% Created in   : 2012-2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update  : 01/07/2014 % Update this date manually
% Last used    : 01/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    return;
end

if nargin < 2
    level = 0;
end

fn = fieldnames(thestruct);

if length(fn) ~= 0
    disp(sprintf('----------------------'))
    disp(sprintf('type\tsize\tfield name'))
    disp(sprintf('----------------------'))
end

for n = 1:length(fn)
    % tabs = '';
    class_var   = eval(['class(thestruct.' fn{n} ');']);
    [a,b]       = eval(['size(thestruct.' fn{n} ');']);
    tabs = sprintf([class_var '\t%0.f x %0.f\t'],a,b);
    for m = 1:level
        tabs = [tabs '    '];
    end
    disp([tabs fn{n}])
    fn2 = getfield(thestruct,fn{n});
    if isstruct(fn2)
        showstruct(fn2, level+1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
