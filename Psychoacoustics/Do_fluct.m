function outs = Do_fluct(insig,fs,dBFS)
% function outs = Do_fluct(insig,fs,dBFS)
%
% 1. Description:
%       Fluctuation implemented by Chalupper and provided to PsySound team.
%       Audio files are calibrated to have 100 dB SPL for 0 dBFS rms. 
%
% 2. Stand-alone example:
%       Do_fluct;
% 
% 3. Additional info:
%       Tested cross-platform: No
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 22/08/2014
% Last update on: 16/02/2015 
% Last use on   : 16/02/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    dBFS = 90;
end

if nargin == 0
    [f1 f2] = uigetfile(Get_TUe_paths('outputs'));
    filename= [f2 f1];
    [x fs]  = Wavread(filename);
    insig   = x;
    dBFS    = input('Enter dB SPL value corresponding to 0 dBFS: ');
end

lvl = dBFS - 114; % To corect from Fastl's calibration to Chalupper's calibration
insig = From_dB(lvl)*insig;
fprintf('%s.m: signal level to be considered approximately as %.1f dB SPL\n',mfilename,rmsdb(insig)+114);

if fs ~= 44100
    error('Model not validated')
end

[N, main_N, spec_N]=ch_dlm(insig);

%calculate loudness fluctuation
[lf lediff outs] = ch_fluct(main_N);

% FS = ( 0.008*sum(outs.le_max-outs.le_min) )/( fmod/4 + 4/fmod )

if nargout == 0
    figure;
    plot(1:24, outs.le_max), hold on
    plot(1:24, outs.le_min), grid on
    ylabel('L_E [dB]')
    xlabel('Critical-band rate [Bark]')
    
    fmod = input('Give fmod value: ');
    FS = ( 0.008*sum(outs.le_max-outs.le_min) )/( fmod/4 + 4/fmod )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
