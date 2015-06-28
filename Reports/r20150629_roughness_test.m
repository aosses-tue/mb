function r20150629_roughness_test
% function r20150629_roughness_test
%
% 1. Description:
%
% 2. Stand-alone example:
%       r20150629_roughness_test;
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 28/06/2015
% Last update on: 28/06/2015 % Update this date manually
% Last use on   : 28/06/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
bDiary = 0;
Diary(mfilename,bDiary);

dir_files = [Get_TUe_paths('outputs') 'Fastl2007_test_20141126' delim]; % These files are calibrated at 90 dB
files = {   'rough_ref.wav', ... % 1
            'rough_test_fc_1000_AM_m_100_fmod_030Hz.wav', ...
            'rough_test_fc_1000_AM_m_100_fmod_050Hz.wav', ...
            'rough_test_fc_1000_AM_m_100_fmod_100Hz.wav', ...
            'rough_test_fc_1000_AM_m_100_fmod_150Hz.wav'; ...
            'rough_test_fc_0125_AM_m_100_fmod_070Hz.wav', ... % 'reference'
            'rough_test_fc_0125_AM_m_100_fmod_030Hz.wav', ...
            'rough_test_fc_0125_AM_m_100_fmod_050Hz.wav', ...
            'rough_test_fc_0125_AM_m_100_fmod_100Hz.wav', ...
            'rough_test_fc_0125_AM_m_100_fmod_150Hz.wav'; ...
            'rough_test_fc_8000_AM_m_100_fmod_070Hz.wav', ... % 'reference'
            'rough_test_fc_8000_AM_m_100_fmod_030Hz.wav', ...
            'rough_test_fc_8000_AM_m_100_fmod_050Hz.wav', ...
            'rough_test_fc_8000_AM_m_100_fmod_100Hz.wav', ...
            'rough_test_fc_8000_AM_m_100_fmod_150Hz.wav'};

N = 8192;
fc = [1000 125 8000];
for i = 1:3

    [xref fs]= Wavread([dir_files files{i,1}]); xref = From_dB(-10)*xref;
    [x1 fs] = Wavread([dir_files files{i,2}]); x1 = From_dB(-10)*x1;
    [x2 fs] = Wavread([dir_files files{i,3}]); x2 = From_dB(-10)*x2;
    [x3 fs] = Wavread([dir_files files{i,4}]); x3 = From_dB(-10)*x3;
    [x4 fs] = Wavread([dir_files files{i,5}]); x4 = From_dB(-10)*x4;

    Rref    = Roughness_offline(xref(1:N),fs);
    R(i,1)    = Roughness_offline(x1(1:N),fs);
    R(i,2)    = Roughness_offline(x2(1:N),fs);
    R(i,3)    = Rref;
    R(i,4)    = Roughness_offline(x3(1:N),fs);
    R(i,5)    = Roughness_offline(x4(1:N),fs);

    RDref   = Roughness_Duisters_offline(xref(1:N),fs);
    if i == 1
        R0 = RDref;
    end
    RD(i,1)   = Roughness_Duisters_offline(x1(1:N),fs);
    RD(i,2)   = Roughness_Duisters_offline(x2(1:N),fs);
    RD(i,3)   = RDref;
    RD(i,4)   = Roughness_Duisters_offline(x3(1:N),fs);
    RD(i,5)   = Roughness_Duisters_offline(x4(1:N),fs);

    RDref2  = Roughness_Duisters_offline(xref(1:2*N),fs,2*N);
    if i == 1
        R02 = RDref2;
    end
    RD2(i,1)  = Roughness_Duisters_offline(x1(1:2*N),fs,2*N);
    RD2(i,2)  = Roughness_Duisters_offline(x2(1:2*N),fs,2*N);
    RD2(i,3)  = RDref2;
    RD2(i,4)  = Roughness_Duisters_offline(x3(1:2*N),fs,2*N);
    RD2(i,5)  = Roughness_Duisters_offline(x4(1:2*N),fs,2*N);

    fmodtest = [30 50 70 100 150];
    figure; % Daniel1997, Fig 3.c
    plot(fmodtest,R(i,:), fmodtest,RD(i,:)/R0, fmodtest,RD2(i,:)/R02);
    legend('Daniel','Duisters','Duisters longer')
    title(['f_c=' num2str(fc(i)) ' [Hz]'])
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
