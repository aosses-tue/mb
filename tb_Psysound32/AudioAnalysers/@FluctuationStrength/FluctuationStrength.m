function obj = FluctuationStrength(varargin)
% function obj = FluctuationStrength(varargin)
%
% 1. Description:
%       FLUCTUATIONSTRENGTH Constructor
%
% 2. Additional info:
%       Tested cross-platform: Yes
%
% Author        : PsySound Team
% Downloaded in : 2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Last update in: 2014 % Update this date manually
% Last use on   : 10/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = struct;

switch nargin
 case 0
  % Default Constructor
  % Inherit from the Analyser base class
  base = Analyser();

  obj = class(obj, 'FluctuationStrength', base);

 case 1
  % Copy Constructor
  % if single argument of class RoughnessDW, return it
  arg1 = varargin{1};
  if isa(arg1, 'FluctuationStrength')
    obj = arg1;
  elseif isstruct(arg1)
    % This should be a file handle
    base = Analyser(arg1);
    
    obj = class(obj, 'FluctuationStrength', base);
  
    % Set window length
    fs  = get(obj, 'fs');
    wl  = fluctuation(fs, 1);
    obj = set(obj, 'windowLength', wl);
    
  else
    error('Fluctuation Strength: Invalid Argument type');
  end
  
 otherwise
  error('Fluctuation Strength: Invalid number of input arguments')
end

% Set name
obj = set(obj, 'Name', 'Fluctuation Strength');

% Set window function
obj = set(obj, 'windowFunc', 'Blackman');

% Specify analyser type
obj = set(obj, 'type', 'Psychoacoustic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end FLUCTUATIONSTRENGTH constructor
