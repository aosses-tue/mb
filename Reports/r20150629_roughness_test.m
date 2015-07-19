function outs = r20150629_roughness_test;
% function outs = r20150629_roughness_test;
%
% 1. Description:
%       It runs the a validation using roughness stimuli.
%       In case the output outs is requested, the validation is not run but
%       some metadata is returned: filenames
% 
% 2. Stand-alone example:
%       r20150629_roughness_test; % it runs the validation
% 
%       outs = r20150629_roughness_test; % it gets the filenames
% 
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 28/06/2015
% Last update on: 29/06/2015 % Update this date manually
% Last use on   : 29/06/2015 % Update this date manually
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
if nargout == 1
    outs.filenames  = files;
    outs.dir_where  = dir_files; 
    return;
end

N = 8192;
fc = [1000 125 8000];
                        % 1 2  3 4
% Duisters_opts.switches = [1 1 2 3];
Duisters_opts.switches = [1 1 12 3];
% Duisters_opts.switches = [2 1 12 3]; THIS IS THE BEST SO FAR
% Duisters_opts.switches = [2 1 2 3]; % Problems at low and high frequencies
% Duisters_opts.switches = [2 1 10 3]; % Worst at low but better at high frequencies

for i = 1:3

    [xref fs]= Wavread([dir_files files{i,1}]); xref = From_dB(-10)*xref;
    [x1 fs] = Wavread([dir_files files{i,2}]); x1 = From_dB(-10)*x1;
    [x2 fs] = Wavread([dir_files files{i,3}]); x2 = From_dB(-10)*x2;
    [x3 fs] = Wavread([dir_files files{i,4}]); x3 = From_dB(-10)*x3;
    [x4 fs] = Wavread([dir_files files{i,5}]); x4 = From_dB(-10)*x4;

    xref = resample(xref,48000,fs);
    x1 = resample(x1,48000,fs);
    x2 = resample(x2,48000,fs);
    x3 = resample(x3,48000,fs);
    x4 = resample(x4,48000,fs);
    fs = 48000;
    
    Rref      = Roughness_offline(xref(1:N),fs);
    R(i,1)    = Roughness_offline(x1(1:N),fs);
    R(i,2)    = Roughness_offline(x2(1:N),fs);
    R(i,3)    = Rref;
    R(i,4)    = Roughness_offline(x3(1:N),fs);
    R(i,5)    = Roughness_offline(x4(1:N),fs);

    idx = 3;
    RDref   = Roughness_Duisters_offline(xref(1:2*N),fs,N,Duisters_opts);
    if i == 1
        try
            R0 = RDref(idx);
            R0first = RDref(1);
        catch
            R0 = RDref(1);
        end
    end
    RD_tmp  = Roughness_Duisters_offline(x1(1:2*N),fs,N,Duisters_opts);
    RD(i,1)   = RD_tmp(idx); RDfirst(i,1) = RD_tmp(1);
    
    RD_tmp    = Roughness_Duisters_offline(x2(1:2*N),fs,N,Duisters_opts);
    RD(i,2)   = RD_tmp(idx); RDfirst(i,2) = RD_tmp(1);
    
    RD(i,3)   = RDref(idx); RDfirst(i,3) = RDref(1);
    
    RD_tmp    = Roughness_Duisters_offline(x3(1:2*N),fs,N,Duisters_opts);
    RD(i,4)   = RD_tmp(idx); RDfirst(i,4) = RD_tmp(1);
    
    RD_tmp    = Roughness_Duisters_offline(x4(1:2*N),fs,N,Duisters_opts);
    RD(i,5)   = RD_tmp(idx); RDfirst(i,5) = RD_tmp(1);

    fmodtest = [30 50 70 100 150];
    figure; % Daniel1997, Fig 3.c
    plot(fmodtest,R(i,:), fmodtest,RD(i,:)/R0); %, fmodtest,RD2(i,:)/R02);
    legend('Daniel','Duisters')%,'Duisters longer')
    title(['f_c=' num2str(fc(i)) ' [Hz]'])
end
RD = RD / R0;
RDfirst = RDfirst / R0first;

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
