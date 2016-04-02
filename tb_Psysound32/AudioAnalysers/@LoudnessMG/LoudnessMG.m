function obj = LoudnessMG(varargin)
% function obj = LoudnessMG(varargin)
% 
% 1. Description:
%       This is the constructor for LoudnessMG.
%       This code implements many of the features of PsySound2, but is suitable
%       ONLY for a sampling frequency fs = 44100 Hz. 
%      This code has been adapted from the original PsySound2 Pascal.
% 
%      The implementation includes the use of the compact spectrum to provide 
%      a computationally efficient solution to loudness calculation.
%      The outputs include:
%        * Loudness (Moore et al static model)
%        * Time varying specific loudness pattern,
%        * Time averaged specific loudness pattern,
%        * Sharpness (Zwicker & Fastl and Aures versions, using the specific 
%          loudness pattern from Moore et al.), 
%        * Volume (vol)
%        * Spectral dissonance (Hutchinson & Knopoff and Sethares),
%        * Tonal dissonance (H&K and Sethares)
%        
% Code by Densil Cabrera and Sam Ferguson
% Some additional comments by Alejandro Osses V., TU/e Eindhoven 2014-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = struct('levelOffset', 0);

switch nargin
 case 0
  % Default Constructor
  % Inherit from the Analyser base class
  base = Analyser();
  
  obj = class(obj, 'LoudnessMG', base);

 case 1
  % Copy Constructor
  % if single argument of class LoudnessMG, return it
  arg1 = varargin{1};
  if isa(arg1, 'LoudnessMG')
    obj = arg1;
  elseif isstruct(arg1)
    % This should be a file handle
    base = Analyser(arg1);
    
    obj = class(obj, 'LoudnessMG', base);
  
  else
    error('LoudnessMG: Invalid Argument type');
  end
  
 otherwise
  error('LoudnessMG: Invalid number of input arguments')
end

% Set name
obj = set(obj, 'Name', 'Loudness (MG & B PsySound2)');

% Set default Overlap, Window size and Windowing function
ov.size = 75;
ov.type = 'percent';

obj = set(obj, 'overlap', ov);
obj = set(obj, 'windowLength', 4096);
obj = set(obj, 'windowFunc', 'Blackman');
% Note - PsySound2 used a hanning window function. However, Blackman is
% probably better.

% Specify analyser type
obj = set(obj, 'type', 'Psychoacoustic');

obj = setlevelOffset(obj, 104.11);
% obj = setlevelOffset(obj, 92.225);

end % LoudnessMG Constructor
