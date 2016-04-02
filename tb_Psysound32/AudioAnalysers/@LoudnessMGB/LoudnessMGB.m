function obj = LoudnessMGB(varargin)
% function obj = LoudnessMGB(varargin)
%
% 1. Description:
%       This is the constructor for LoudnessMGB: 'Model for estimating the 
%       loudness of time-varying sounds' following Moore, Glasberg and Baer. 
%     
%       'THIS ANALYSER IS SLOW, therefore test a short duration recording first.
% 
%     References:
%       B.R. Glasberg and B.C.J. Moore. 1990. 'Derivation of Auditory Filter 
%           Shapes from Notched Noise Data. Hearing Research, 47: 103-137.
% 
%       B.C.J. Moore, B.R. Glasberg and T. Baer. 1997. A Model for the 
%           Prediction of Thresholds, 'Loudness, and Partial Loudness. 
%           Journal of the Audio Engineering Society, 45(4): 224-240.
%
%       B.R. Glasberg and B.C.J. Moore. 2002. A Model of Loudness Applicable 
%           to Time-Varying Sounds. Journal of the Audio Engineering 
%           Society, 50(5): 331-342. 
% 
%       B.R. Glasberg and B.C.J. Moore. 2006. Prediction of absolute 
%           thresholds and equal-loudness contours using a modified loudness 
%           model. J. Acoust. Soc. Am. 120: 585-588.
% 
%       B.C.J. Moore and B.R. Glasberg. 2007. Modelling Binaural Loudness
%           J. Acoust. Soc. Am. 121: 1604-1612
%
% Code by D Cabrera, D Lee & S Ferguson.
% Some additional comments by Alejandro Osses V., TU/e Eindhoven 2014-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = struct('filterMethod',1);

switch nargin
 case 0
  % Default Constructor
  % Inherit from the Analyser base class
  base = Analyser();

  obj = class(obj, 'LoudnessMGB', base);

 case 1
  % Copy Constructor
  % if single argument of class LoudnessMGB, return it
  arg1 = varargin{1};
  if isa(arg1, 'LoudnessMGB')
    obj = arg1;
  elseif isstruct(arg1)
    % This should be a file handle
    base = Analyser(arg1);
    
    obj = class(obj, 'LoudnessMGB', base);
  
  else
    error('LoudnessMGB: Invalid Argument type');
  end
  
 otherwise
  error('LoudnessMGB: Invalid number of input arguments')
end

% Specify analyser type
obj = set(obj, 'type', 'Raw');

% Set default windowlength
numSamples = get(obj,'samples') + mod(get(obj,'samples'),2);
obj = set(obj, 'windowLength', numSamples);

% set output rate to 100 - if synchronise is used this will be changed
% This considerably decreases storage space necessary.
% obj = set(obj, 'outputDataRate', 100);

% Set stereo mode
obj = set(obj, 'multiChannelSupport', true);

% Set name
obj = set(obj, 'Name', 'Loudness (Moore, Glasberg and Baer)');

% end LoudnessMGB constructor
