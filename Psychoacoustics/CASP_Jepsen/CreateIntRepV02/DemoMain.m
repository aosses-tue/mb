%--------------------------------------------------------------------------
% DemoMain.m
%--------------------------------------------------------------------------
% 1. Description:
%   runs preprocessing stages of pemo resp. casp and calculates and plots
%   the internal representation of an arbitrary input signal
%
% usage
%   DemoMain(model, pmode, wav)
%
% input
%   model       : string specifying if pemo or casp is used
%   pmode       : string specifying what is plotted, possible options are:
%                   - input         : input signal
%                   - basilar       : output basilar-membrane filterbank
%                   - basilarrms    : rms(x(t)) at output b-m filterbank
%                   - hc            : output after hair cell transformation
%                                     stage
%                   - adapt         : output adaptation loops
%                   - IntRep        : internal representation of input,
%                     includes several options depending on the range of
%                     basilar-membrane and modulation filterbank, please
%                     browse PlotIntRep.m for details
%  wav          : string containing the filename of the input signal in
%                 wav-format
%
% 2. Examples:
%       model = 'casp';
%       DemoMain(model);
% 
% version 1.0
%   20/01/2013, C.T. Iben
% version 2.0
%   29/01/2013, C.T. Iben

function DemoMain(model, pmode, wav)

%% if no input arguments are set the demo version will be run
if nargin < 3, wav   = 'x.wav'; end % test stimulus1: 'x.wav', test stimulus2: 'S10M_duhd.wav'
if nargin < 2, pmode = 'all'; end
if nargin < 1, model = 'pemo'; end

fprintf('Tool is running with %s in plotmodus %s with signal %s\n',model,pmode,wav)

%% set matlab path
startupDemoMain; % name changed by AO startup was shadowed by my startup file

%% load input signal
[x, fs] = Wavread(wav);

%% runs preprocessing of the specified model and calculates the internal representation
switch model
    case 'pemo'
        [IntRep BM MF] = PemoPreProc(x, fs);
    case 'casp'
        [IntRep BM MF] = CaspPreProc(x, fs);
    otherwise
        error('illegal model modus')
end

%% plots output of pre-processing stages and/or internal representation of the input signal please browse PlotIntRepMain.m for details
PlotIntRepMain(IntRep, BM, MF, fs, pmode)

%% assigns structure IntRep to workspace
assignin('base','IntRep',IntRep);
assignin('base','BM',BM);
assignin('base','MF',MF);
disp('Model structures IntRep, BM and MF now in current workspace available')

% eof
