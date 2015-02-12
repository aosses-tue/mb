function l0P392_20150209
% function l0P392_20150209
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 08/02/2015
% Last update on: 08/02/2015 % Update this date manually
% Last use on   : 08/02/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

% dir_scripts = ''; % 'here' in case is empty
% dir_wav     = '~/Documenten/Documenten-TUe/09-Training+activities/2015-Q3-Advanced-perception/Assigment1_Sound_files/';

% inputfile - 't54.wav'
lowpassfilter;

% inputfile - 't54.wav'
quant;

%% inputfile = 'trumpet.wav';
% Quantisation to 10 bits:
quant_example1;

% Quantisation per segment (1-sec length) between 10 and 12 bits
quant_example2;

% Quantisation per segment (50-ms length): looks for No of bits producing 
% an SNR of approx. 30 dB
quant_example3;

insig1 = 'Track-2-L-excerpt-18.wav';

outsig1 = [Delete_extension(insig1,'wav') '-quant-SNR-30-dB.wav'];
quant_example3_s(insig1,outsig1,30);

outsig1 = [Delete_extension(insig1,'wav') '-quant-SNR-15-dB.wav'];
quant_example3_s(insig1,outsig1,15);

resample;

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
