function l0P392_20150225
% function l0P392_20150225
%
% 1. Description:
%	For lecture of the advanced perception course (Speech perception 1)
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 25/02/2015
% Last update on: 25/02/2015 % Update this date manually
% Last use on   : 25/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

dir = '~/Documenten/Documenten-TUe/09-Training+activities/2015-Q3-Advanced-perception/College/';
dirmedia    = [dir 'media' delim];
dirout      = [dir 'output' delim];
Mkdir(dirout);

bLPC = 1;
bFigures4lesson = 1;

if bFigures4lesson
    [x fs] = Wavread([dirmedia 't49.wav']);
    t = ( 0:length(x)-1 )/fs;
    figure;
    plot(t,x), grid on;
    title(sprintf('''Do you want to create a new document'', fs=%.0f [Hz]',fs))
    xlim(minmax(t));
    xlabel('t (s)')
    ylabel('Amplitude')

    opt.format = 'png';
    Saveas(gcf,[dirout 't49-waveform'],opt);
end

if bLPC == 1
    fname1 = [dir 'rl001/rl001-car.wav'];
    fname2 = [dir 'rl001/rl001-park.wav'];
    fname3 = [dir 'sb001/sb001-car.wav'];
    fname4 = [dir 'sb001/sb001-park.wav'];

    [x1 fs]  = Wavread(fname1);
    [x2 fs]  = Wavread(fname2);
    [x3 fs]  = Wavread(fname3);
    [x4 fs]  = Wavread(fname4);

    N = 8192;
    Get_LPC(x1,fs,N);
    % Formant 1 Frequency 674.5
    % Formant 2 Frequency 1296.0
    % Formant 3 Frequency 2246.9
    % Formant 4 Frequency 3483.4
    % Formant 5 Frequency 4222.4

    Get_LPC(x2,fs,N);
    % Formant 1 Frequency 664.7
    % Formant 2 Frequency 1135.7
    % Formant 3 Frequency 2231.3
    % Formant 4 Frequency 3255.0
    % Formant 5 Frequency 4205.5

    Get_LPC(x3,fs,N);
    % Formant 1 Frequency 707.0
    % Formant 2 Frequency 999.5
    % Formant 3 Frequency 1765.2
    % Formant 4 Frequency 3392.0
    % Formant 5 Frequency 4432.3

    Get_LPC(x4,fs,N);
    % Formant 1 Frequency 549.7
    % Formant 2 Frequency 924.4
    % Formant 3 Frequency 1455.9
    % Formant 4 Frequency 2617.3
    % Formant 5 Frequency 3452.8
    % Formant 6 Frequency 4434.2
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
