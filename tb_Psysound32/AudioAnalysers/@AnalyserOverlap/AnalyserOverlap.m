function obj = AnalyserOverlap(varargin)
% function obj = AnalyserOverlap(varargin)
%
% 1. Description:
% 
%   ANALYSER Constructor.  This is the base class for all the Analysers
% 
%   OBJ = Analyser          Creates an empty Analyser
%   OBJ = Analyser(Obj)     Returns ObjIn
%   OBJ = Analyser(Obj, fH) Creates a new Analyser object and
%                           initialises its fields using the file
%                           handle, fH
%   Enhancement: Consider datatyping the fields to enforce stronger type 
%   checking.
% 
% 2. Additional info:
%       Tested cross-platform: Yes
%
% Author        : PsySound Team
% Downloaded in : 2014
% Modified by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Adapted from  : Analyser.m
% Last update in: 10/09/2014 % Update this date manually
% Last use on   : 10/09/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
 case 0
  % if no input arguments, create a default object
  obj.Name         = 'Base Analyser';
  obj.filename     = '';
  obj.fs           = 0;
  obj.bits         = 0;
  obj.samples      = 0;
  obj.channels     = 0;
  obj.windowLength = -1;
  obj.windowFunc   = 'rect';
  obj.overlap      = struct('size', 4096, 'type', 'samples');
  obj.output       = {};
  obj.SummaryOutput = {}; %Added New release
  obj.OptionStr = {}; %Added New release

  % 
  % This field specifies the Analyser type in terms of how the data
  % is processed. i.e. in windows like FFT and psychoacoustics
  % algorithms or via time domain filtering eg. SLM. The latter has
  % the property that the output is a timeseries with the same
  % length as the data.  Also affects synchronisation.
  %
  obj.type = '';
  
  % This allows for samples from the beginning & end of each
  % window to be thrown away.  See Hilbert for an example. The
  % default, of course, is zero.
  obj.windowOffset = 0;
  
  % Rate of output data
  obj.outputDataRate = 0;
  
  % Number of output samples
  obj.outputSamples = 0;
  
  % Signify synchronisation
  obj.synch = false;
 
  % Specifies multichannel support
  obj.multiChannelSupport = false;
  
  % Specify class name
  obj = class(obj, 'AnalyserOverlap');
 
 case 1
  % if single argument of class Analyser, return it
  % if single argument of struct, then initialise fields using the
  %                                    file handle
  arg1 = varargin{1};
  if isa(arg1, 'AnalyserOverlap') 
    obj = arg1;
  elseif isstruct(arg1) 
    fh  = arg1;
    obj = AnalyserOverlap;

    % Call nested function
    initParams();
  else
    error('AnalyserOverlap: Invalid argument type')
  end 
  
 otherwise
  error('AnalyserOverlap: Invalid number of input arguments')
end

  function initParams
  % xxx - This could be streamlined somewhat by using a FileAttribute
  %      object
    obj = set(obj, 'filename', fh.name);
    obj = set(obj, 'fs',       fh.sampleRate);
    obj = set(obj, 'bits',     fh.bitsPerSample);
    obj = set(obj, 'channels', fh.channels);
    obj = set(obj, 'samples',  fh.samples);
  end % initParams
end % Analyser constructor

% Field descriptions

% filename
% fs
% samples
% bits
% channels
%
% These fields are taken from the fileHandle directly. They're saved here
% to avoid confusion resulting from trying to associate fileHandles with
% results after the analysis process. 

% windowlength
% overlap
%
% These fields are sometimes user selectable, depending on the amount of
% resolution that people require, and sometimes they are set by the
% algorithm.

% output
%
% This is a cell array used to save the dataObjects in when the process
% command is called. 
