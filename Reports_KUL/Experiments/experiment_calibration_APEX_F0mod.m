function apex=experiment_calibration_APEX_F0mod(p, targetpath, filename)
% function apex=experiment_calibration_APEX_F0mod(p, targetpath, filename)
%
% Returns calibration signal and apex experiment XML code 
%
% Adapted from score_calibration by Alejandro Osses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filename = wav-file name of calibration tone

p = Ensure_field(p,'speechmaterial_uri',[]);

switch p.speechmaterial
    case 'lilliput'
        profile='spin-lilliput_ltass';
        % calibrationlevel=60;
        % siglevel=64;
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'weightednoise_nosilence28.4.wav'];
    case 'LISTf'
        calibrationlevel=60;
        siglevel=64;
        profile='spin-LISTvrouw_ltass';
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'noise/wivineruis.wav'];
    case 'LISTm'
        calibrationlevel=60;
        siglevel=64;
        profile='spin-LISTman_ltass';
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'noise/jwruis.wav'];
    case 'VlMatrix'
        calibrationlevel=60;
        siglevel=50;
        profile='VlMatrix_ltass'; 
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'VlMatrixnoise_ltass.wav'];
    case 'EsMatrix'
        calibrationlevel=60;
        siglevel=50;
        profile='VlMatrix_ltass'; 
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'VlMatrixnoise_ltass.wav'];    
    otherwise
        warning('Default speech material...')
        calibrationlevel=90; % play calibration signal at 90dBSPL
        siglevel=60;
        profile=sprintf('score_%d', p.calibration);
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'sin_1kHz_' num2str(siglevel-p.calibration) 'dBFS.wav'];
end

fs  = 44100;
f   = 1000;
len = 5;
signal = sin(2*pi*f*[0:len*fs-1]/fs);
% signal = signal/rms(signal) * 10^((siglevel-q.calibration+transcorr)/20);

dbid='db_calib';
stimid='stim_calib';

pc.calibration_amplitude=calibrationlevel;
apex.calibration    = a3calibration(profile, stimid, cal_stim_id, siglevel, pc);
apex.datablock      = a3datablock(dbid, filename, 'soundcard');
apex.stimulus       = a3stimulus(stimid, dbid);

apex.stimid=stimid;
apex.dbid=dbid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end