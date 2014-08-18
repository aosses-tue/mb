function Generate_reference_sounds
% function Generate_reference_sounds
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: Yes
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 16/08/2014 % Update this date manually
% Last use on   : 16/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bZeroPad = 1;
bSave = 1;

% Common params
fc  = 1000;
T   = 1/fc;
fs  = 44100;
dur = 4;
sig = Create_sin(fc,dur,fs,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Loudness
lvl = 40;
lvlAMT  = lvl + 10; % Zwicker's correction
y = setdbspl(sig,lvlAMT);

if bZeroPad
    yzero = Zero_padding(y,200e-3,fs);
end

if bSave
    Wavwrite(yzero,fs,[Get_TUe_paths('outputs') 'ref_loud']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Sharpness % 
fmod    = 70;
m       = 100;
option  = 'm';
lvl     = 70;

sigbp = Wavread('~/Documenten/MATLAB/outputs/track_04_t09_bp.wav');

start_phase = pi/2; % begin in maximum. Use -pi/2 to begin in minimum
% ch_am(sig,fmod,m,option,fs,start_phase); % uncomment this to plot
y = ch_am(sigbp,fmod,m,option,fs,start_phase);

lvlAMT  = lvl + 10; % Zwicker's correction
y = setdbspl(y,lvlAMT);

if bZeroPad
    yzero = Zero_padding(y,200e-3,fs);
end

if bSave
    Wavwrite(yzero,fs,[Get_TUe_paths('outputs') 'ref_sharp']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Fluctuation strength
fmod    = 4;
m       = 100;
option = 'm';
lvl     = 60;

start_phase = pi/2;
y = ch_am(sig,fmod,m,option,fs,start_phase);

lvlAMT  = lvl + 10; % Zwicker's correction
y = setdbspl(y,lvlAMT);

if bZeroPad
    yzero = Zero_padding(y,200e-3,fs);
end

if bSave
    Wavwrite(yzero,fs,[Get_TUe_paths('outputs') 'ref_fluct']);
end

% 4. Roughness
fmod    = 70;
m       = 100;
option = 'm';
lvl     = 60;
start_phase = pi/2;
y = ch_am(sig,fmod,m,option,fs,start_phase);

lvlAMT  = lvl + 10; % Zwicker's correction
y = setdbspl(y,lvlAMT);

if bZeroPad
    yzero = Zero_padding(y,200e-3,fs);
end

if bSave
    Wavwrite(yzero,fs,[Get_TUe_paths('outputs') 'ref_rough']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
