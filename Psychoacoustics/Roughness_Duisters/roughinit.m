function [] = roughinit(fs, switches);

% [] = roughinit(fs, switches);
%
% Initialises the roughness model by creating variables and saving them
% to the file roughinit.mat.
% Variables are created depending on sample frequency ’fs’, and the
% selected module implementations defined in ’switches’.

global numH7 denH7 numH14 denH14 numH30 denH30 numH36 denH36 numH66 denH66

if nargin < 2
    switches = [2213];
end

if nargin < 1
    fs = 48000;
end

t = 0:1/fs:0.6;
if switches(4) == 3
    disp('calculating filter weights...');
    [numH7, denH7, numH14, denH14, numH30, denH30, numH36, denH36, numH66, denH66] ...
    = fir_hweight(fs);
    twin = 0.5;
    dN = [round(0.3*fs) round(0.3*fs+fs/70)];
    ref1 = (1 + sin(2 * pi * 70 * t)) .* sin(2 * pi * 1e3 * t);
    win = sin(0.5 * pi * 50 .* (0:1/fs:20e-3)).^2;
    ref1 = [ref1(1:length(win)) .* win ref1(length(win)+1:length(ref1)-length(win)) ...
    ref1(length(ref1)-length(win)+1:length(ref1)) .* fliplr(win)];
    disp(’calculating asper reference...’);
    R01 = roughcalc(ref1, 60, fs, twin, dN, 1.6, switches);
    save roughinit fs numH7 denH7 numH14 denH14 numH30 denH30 numH36 denH36 numH66 denH66 switches R01
else
    twin = 0.6;
    dN = [floor(0.3*fs) floor(0.5*fs)];
    ref1 = (1 + sin(2 * pi * 70 * t)) .* sin(2 * pi * 1e3 * t);
    win = sin(0.5 * pi * 50 .* (0:1/fs:20e-3)).^2;
    ref1 = [ref1(1:length(win)) .* win ref1(length(win)+1:length(ref1)-length(win)) ...
    ref1(length(ref1)-length(win)+1:length(ref1)) .* fliplr(win)];
    disp('calculating asper reference...');
    R01 = roughcalc(ref1, 60, fs, twin, dN, 1.6, switches);
    save roughinit fs switches R01
end
