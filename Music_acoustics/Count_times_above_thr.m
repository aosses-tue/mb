function [Counter params] = Count_times_above_thr(filename)
% function [Counter params] = Count_times_above_thr(filename)
%
% 1. Description:
%   Measures rotation times related to the voice of the dragon instrument. 
%   Period 'i' starts at params.tCounter(i), given in seconds.
% 
%       params.ti   - initial time for each period [s]
%       params.T    - each period duration [s]
%         
% 2. Additional info:
%   Tested cross-platform: No
%
% 3. Stand-alone example:
%   Count_times_above_thr;
% 
% Programmed by Alejandro Osses, HIT, TU/e, the Netherlands, 2014
% Created on : 26/05/2014
% Last update: 25/06/2014 % Update this date manually
% Last used  : 23/07/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    path.db_voice_of_dragon = Get_TUe_paths('db_voice_of_dragon');
    [f1 f2] = uigetfile([path.db_voice_of_dragon delim '02-Wav-files' delim '1 referentie' delim '*3.wav'],'*.wav');
    filename = [f2 f1];
end

thr_up      = 0.5;
thr_down    = -0.1;
bCount      = 1;
Counter     = 0;
tCounter    = 0;
waittime    = 0.15;
[x fs] = wavread(filename);
t = (0:length(x)-1)/fs;
for i=1:length(x)
    if bCount == 1 & x(i)>= thr_up
        Counter = Counter + 1;
        tCounter = [tCounter; t(i)];
        bCount = 0;
    elseif (t(i)-tCounter(end))>=waittime
        if x(i)<= thr_down
            bCount = 1;
        end
    end
end

tCounter(1)     = [];

params.T        = diff(tCounter);

tCounter(end)   = []; % last period is incomplete

params.ti       = tCounter;
params.N        = length(params.T); % should be the same than 'Counter'

if nargin==0
    figure;
    plot(t,x)
    title([num2str(Counter) ' rotations found. Ti = ' num2str(tCounter(1))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end