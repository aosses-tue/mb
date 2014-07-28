function Get_filtered_signals(dir_in, dir_out)
% function Get_filtered_signals(dir_in, dir_out)
%
% 1. Description:
%
% 2. Additional info:
%   Tested cross-platform: No
%
% 3.1 Stand-alone example:
%       dir_in  = [Get_TUe_paths('db_voice_of_dragon') '02-Wav-files\04-Wav-files-calibrated-44.1kHz\'];
%       dir_out = [Get_TUe_paths('db_voice_of_dragon') '02-Wav-files\05-Wav-files-calibrated-44.1kHz-filtered\'];
%       Get_filtered_signals(dir_in, dir_out);
% 3.2
%       dir_in  = [Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted\04-Wav-files-calibrated-44.1kHz\'];
%       dir_out = [Get_TUe_paths('db_voice_of_dragon') '03-Wav-files-predicted\05-Wav-files-calibrated-44.1kHz-filtered\'];
%       Get_filtered_signals(dir_in, dir_out);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 3/7/2014
% Last update on: 3/7/2014 % Update this date manually
% Last used on  : 3/7/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    dir_in = uigetdir('Select the directory where the wav-files to be filtered are: ');
    dir_in = [dir_in delim];
end

if nargin < 2
    dir_out = [Get_TUe_paths('outputs') 'new_wav' delim];
end

Mkdir(dir_out);

Diary(mfilename,dir_out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extra.bExtension = 0; % To delete extension
filenames        = Get_filenames(dir_in,['*.wav'],extra);

BandsPerOctave = 1;
f = Get_OB_freqs(BandsPerOctave);

fcmin = 100;
fcmax = 5000;

fc = f(find(f>fcmin & f<fcmax));
fcut_inf = fc/2; % 1 octave below
fcut_sup = fc*2; % 1 octave above

dBBelow = 20;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Signals are going to be filtered with the following centre freqs, cut-off freqs (inf) and cut-off freqs (sup) [Hz]')
disp(num2str(round(fc)))
disp(num2str(round(fcut_inf)))
disp(num2str(round(fcut_sup)))
disp(['A message is going to be displayed for every filtered signal having an RMS > xrms - ' num2str(dBBelow) ' dB'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

for i = 1:length(filenames)
    disp(filenames{i});
    
    [x fs] = Wavread([dir_in filenames{i}]);
    info.fs = fs;
    xrms = rmsdb(x);
    
    for j = 1:length(fc)
        y = freqfftwhpf( x, info, fcut_inf(j) ); % default fc = 50 Hz;
        y = freqfftwlpf( y ,info, fcut_sup(j) );
        outputfilename = [dir_out filenames{i} '-fc-' num2str(round(fc(j))) '-Hz.wav'];
        yrms = rmsdb(y);
        
        if yrms > xrms - dBBelow
            Wavwrite(y,fs,outputfilename);
        else
            wavwrite(y,fs,outputfilename);
        end
        
    end

end

diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end