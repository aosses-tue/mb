function [y outs] = Generate_test_AMtone(fc,fmod,m_or_d,option,lvl,fs,dur,dur_ramp_ms,start_phase)
% function [y outs] = Generate_test_AMtone(fc,fmod,m_or_d,oprion,lvl,fs,dur,dur_ramp_ms,start_phase)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/03/2015
% Last update on: 05/03/2015 % Update this date manually
% Last use on   : 05/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    fc  = 1000;
end

if nargin < 2
    fmod = 70; % assuming d = 40 dB
end

if nargin < 3
    m_or_d = 40; % assuming d = 40 dB
end

if nargin < 4
    option = 'd'; % 'm'
end

if nargin < 5
    lvl  = 70;
end
lvlAMT  = lvl + 10;

if nargin < 6
    fs  = 44100;
end

if nargin < 7
    dur = 1; % 1 second
end

if nargin < 8
    dur_ramp_ms = 0;
end

if nargin < 9
    start_phase = 0;
end

sig = Create_sin(fc,dur,fs,0);

if      strcmp(option,'d')
    y = ch_am(sig,fmod,m_or_d,option,fs,start_phase);
elseif  strcmp(option,'m')
    y = ch_am(sig,fmod,m_or_d,option,fs,start_phase);
end

y = setdbspl(y,lvlAMT);

ramp2apply = cos_ramp(length(y),fs,dur_ramp_ms);
y    = ramp2apply'.*y;

if      strcmp(option,'d')
    filename = [Get_TUe_paths('outputs') 'test_fc_' Num2str(fc   ,4) '_AM_d_' Num2str(m_or_d,3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
elseif  strcmp(option,'m')
    filename = [Get_TUe_paths('outputs') 'test_fc_' Num2str(fc   ,4) '_AM_m_' Num2str(m_or_d,3) '_fmod_' Num2str(floor(fmod),3) 'Hz_' num2str(lvl) '_dBSPL'];
end

if nargout ~= 0
    Wavwrite(y   ,fs,filename);
else
    t = ( 1:length(y) )/fs;
    figure;
    plot(t,y)
end

outs.filename = filename;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
