function h = Analyse_audio_classic_psychoacoustics(fi1,fi2,fi3,options,stPlot)
% function h = Analyse_audio_classic_psychoacoustics(fi1,fi2,fi3,options,stPlot)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 31/10/2014
% Last update on: 31/10/2014 % Update this date manually
% Last use on   : 31/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    options = [];
end

h = [];

options = Ensure_field(options,'bDoFluct' ,0);
options = Ensure_field(options,'bDoFluct3',0); % Fluctuation for 3 audio signals
options = Ensure_field(options,'bDoPsySound',1);

bDoFluct    = options.bDoFluct;
bDoFluct3   = options.bDoFluct3;

ti          = options.ti;
tf          = options.tf;

%% Load audio data

[x1 fs] = Wavread(fi1);
[x2   ] = Wavread(fi2);
try
    [x3   ] = Wavread(fi3);
    bDoFluct = 0;   
catch
    x3 = [];
end
    
%%
if bDoFluct
    
    dBFS = 100; % Zwicker's calibration

    ttemp = ( 1:length(x1) )/fs;
    idx = find(ttemp > ti-0.5 & ttemp < tf+0.5 );
    x1 = x1(idx);
    x2 = x2(idx);
    
    [h(end+1:end+2) tmp] = Do_fluc_20140930(x1,x2,stPlot,fs);

end

if bDoFluct3
    
    dBFS = 100; % Zwicker's calibration

    ttemp = ( 1:length(x1) )/fs;
    idx = find(ttemp > ti-0.5 & ttemp < tf+0.5 );
    x1 = x1(idx);
    x2 = x2(idx);
    
    h(end+1) = Do_fluc_20141031(x1,x2,x3,stPlot,fs);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
