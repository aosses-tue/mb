function info = Check_Fastl2007(chapter,fig_number)
% function info = Check_Fastl2007(chapter,fig_number)
%
% 1. Description:
%       Tracks ready:
%       38
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       % Chapter 7, figure 8 (Figure 7.8)
%       chapter = 7;
%       fig_number = 8;
%       Check_Fastl2007(chapter,fig_number);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/08/2014
% Last update on: 13/08/2014 % Update this date manually
% Last use on   : 13/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% track_35  -   Figure 10.1

if nargin == 0
    chapter = 10;
end

if nargin < 2
    fig_number = 1;
end

paths = Get_TUe_paths('db_fastl2007');

switch chapter
    case 1
        switch fig_number
            case 1
                filename = 'track_04';
                timesep = [ 1 5 7 11 13 17 ... % 1 kZ, beat, AM tone
                            39 43 45 49 51 55]; % 39 = white; 45 = band-pass, 51 = narrow-band
        end
    case 7 % Just-notiecable sound changes
        switch fig_number
            case 8
                filename = 'track_28';
        end
    case 8
        switch fig_number
            case 12
                filename = 'track_33';
                timesep  = [ 0      0+1 2 2.5 3.5  3.5+1  5.5  6    7, ...
                             9      9+1 10.5 10.8 10.8+1 12.8 13.3 13.6, ...
                            15.6 15.6+1 17.1 17.2 18.2 19.2 19.7 19.8]; % Continue...
        end
    case 9
        switch fig_number
            case 3
                filename = 'track_34';
                timesep = [3 4 7 8 11 14 17 18 21 22]; % not perfect grid
        end
    case 10 % Fluctuation strength
        switch fig_number
            case 1
                filename = 'track_35';
                timesep = [  1  6  8 13 15 20 ...
                            24 29 31 36 38 43 ...
                            47 52 54 59 61 66];
            case 6
                filename = 'track_36';
            case 7
                filename = 'track_37';
        end
    case 11 % Roughness 
        switch fig_number
            case 1
                filename = 'track_38';
                timesep = [1:3:41];
        end
end

[x fs] = Wavread([paths filename]);

t = ( 1:length(x) )/fs;

info.x = x;
info.fs = fs;
info.t = t;
info.timesep = timesep;

figure
plot(t,x);
grid on;
xlabel('Time [s]')
ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
