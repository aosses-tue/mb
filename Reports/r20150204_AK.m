function r20150204_AK
% function r20150204_AK
%
% 1. Description:
%       Little favour that Armin asked me to test the phase cancellation. 
%       Two audio files were generated and presented to Armin presentially.
% 
%       Some concepts:
%           - Random phase
%           - Lateral inhibition (at some point of the conversation)
%           - Gammatone filterbank vs Strube's filterbank. The first was 
%           giving the same filter output for all the phase conditions, while
%           the Strube filterbank was not. That means that there are some
%           nonlinearities that are better addressed by the Strube model.
% 
% 2. Stand-alone example:
%       r20150204;
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 04/02/2015
% Last update on: 04/02/2015 % Update this date manually
% Last use on   : 04/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

output_dir = [Get_TUe_paths('outputs') 'AK' delim];

dur = 4;
fs = 44100;

A = 0.005;
f0 = 6;
n_harm = ceil(2000/f0):floor(8000/f0);
fn = f0*n_harm;

a = -pi;
b = pi;
random_phases = a + (b-a).*rand(length(fn),1);

start_phase = pi/2; % starts in amplitude 1
[y, t]  = Create_sin_phase(f0,start_phase, dur,fs);

yh      = zeros(size(y));
yh_rand = zeros(size(y));

for i = 1:length(fn)
    ytmp  = Create_sin_phase(fn(i),start_phase,dur,fs);
    yh = yh + ytmp;
    ytmp2 = Create_sin_phase(fn(i),start_phase+random_phases(i),dur,fs);
    yh_rand = yh_rand + ytmp2;

end

yt = A*y+A*yh;
yt_rand = A*y+A*yh_rand;

lvl = 70;
yt = setdbspl(yt,lvl);
yt_rand = setdbspl(yt_rand,lvl);

Wavwrite(yt     ,fs,sprintf('%scos-f-phase-0'   ,output_dir));
Wavwrite(yt_rand,fs,sprintf('%scos-f-phase-rand',output_dir));

figure;
plot(t,yt), grid on, hold on
plot(t,yt_rand,'r')


if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
