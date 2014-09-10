function obj = RoughnessDH(varargin)
% function obj = RoughnessDH(varargin)
%
% 1. Description:
%       ROUGHNESSDH Constructor. Overlap is here set to 50%, i.e., 2048 samples
%
% 2. Additional info:
%       Tested cross-platform: Yes
%
% Author        : PsySound Team
% Downloaded in : 2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Adapted from  : RoughnessDW.m
% Last update in: 2014 % Update this date manually
% Last use on   : 10/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = struct;

switch nargin
 case 0
  % Default Constructor
  % Inherit from the Analyser base class
  base = AnalyserOverlap();

  obj = class(obj, 'RoughnessDW', base);

 case 1
  % Copy Constructor
  % if single argument of class RoughnessDW, return it
  arg1 = varargin{1};
  if isa(arg1, 'RoughnessDW')
    obj = arg1;
  elseif isstruct(arg1)
    % This should be a file handle
    base = Analyser(arg1);
    
    obj = class(obj, 'RoughnessDW', base);
  
    % Set window length
    fs  = get(obj, 'fs');
    wl  = roughness(fs, 1);
    obj = set(obj, 'windowLength', wl);
    
  else
    error('RoughnessDW: Invalid Argument type');
  end
  
 otherwise
  error('RoughnessDW: Invalid number of input arguments')
end

% Set name
obj = set(obj, 'Name', 'Roughness (D & W)');

% Set window function
obj = set(obj, 'windowFunc', 'Blackman');

% Specify analyser type
obj = set(obj, 'type', 'Psychoacoustic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end RoughnessDW constructor
