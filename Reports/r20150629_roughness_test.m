function r20150629_roughness_test
% function r20150629_roughness_test
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 28/06/2015
% Last update on: 28/06/2015 % Update this date manually
% Last use on   : 28/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

dir_files = 'D:\Output\Fastl2007_test_20141126\'; % These files are calibrated at 90 dB
files = {   'rough_ref.wav', ... % 1
            'rough_test_fc_1000_AM_m_100_fmod_030Hz.wav', ...
            'rough_test_fc_1000_AM_m_100_fmod_050Hz.wav', ...
            'rough_test_fc_1000_AM_m_100_fmod_100Hz.wav', ...
            'rough_test_fc_1000_AM_m_100_fmod_150Hz.wav'};

N = 8192;
[xref fs]   = Wavread([dir_files files{1}]); xref = From_dB(-10)*xref;
[x1 fs]     = Wavread([dir_files files{2}]); x1 = From_dB(-10)*x1;
[x2 fs]     = Wavread([dir_files files{3}]); x2 = From_dB(-10)*x2;
[x3 fs]     = Wavread([dir_files files{4}]); x3 = From_dB(-10)*x3;
[x4 fs]     = Wavread([dir_files files{5}]); x4 = From_dB(-10)*x4;

Rref = Roughness_offline(xref(1:N),fs);
R(1) = Roughness_offline(x1(1:N),fs);
R(2) = Roughness_offline(x2(1:N),fs);
R(3) = Rref;
R(4) = Roughness_offline(x3(1:N),fs);
R(5) = Roughness_offline(x4(1:N),fs);

RDref = Roughness_Duisters_offline(xref(1:N),fs);
RD(1) = Roughness_Duisters_offline(x1(1:N),fs);
RD(2) = Roughness_Duisters_offline(x2(1:N),fs);
RD(3) = RDref;
RD(4) = Roughness_Duisters_offline(x3(1:N),fs);
RD(5) = Roughness_Duisters_offline(x4(1:N),fs);

RDref2 = Roughness_Duisters_offline(xref(1:2*N),fs,2*N);
RD2(1) = Roughness_Duisters_offline(x1(1:2*N),fs,2*N);
RD2(2) = Roughness_Duisters_offline(x2(1:2*N),fs,2*N);
RD2(3) = RDref2;
RD2(4) = Roughness_Duisters_offline(x3(1:2*N),fs,2*N);
RD2(5) = Roughness_Duisters_offline(x4(1:2*N),fs,2*N);

fmodtest = [30 50 70 100 150];
figure;
plot(fmodtest,R, fmodtest,RD/RDref, fmodtest,RD2/RDref2);
legend('Daniel','Duisters','Duisters longer')

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
