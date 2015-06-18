function apex = experiment_selection_noise_APEX(p)
% function apex = experiment_selection_noise_APEX(p)
%
% Returns calibration signal and apex experiment XML code 
%
% Original file: experiment_calibration_APEX_F0mod.m
% Created on: 18/06/2015
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
        warning('Fix noise as VlMatrix and EsMatrix')
    case 'LISTf'
        calibrationlevel=60;
        siglevel=64;
        profile='spin-LISTvrouw_ltass';
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'noise' delim 'wivineruis.wav'];
        warning('Fix noise as VlMatrix and EsMatrix')
    case 'LISTm'
        calibrationlevel=60;
        siglevel=64;
        profile='spin-LISTman_ltass';
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'noise' delim 'jwruis.wav'];
        warning('Fix noise as VlMatrix and EsMatrix')
    case 'VlMatrix'
        p = Ensure_field(p,'noisefile','VlMatrixnoise_ltass.wav');    
        calibrationlevel=60;
        siglevel=50;
        profile='VlMatrix_ltass'; 
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri p.noisefile];
        warning('Fix noise as VlMatrix and EsMatrix')
    case 'EsMatrix'
        p = Ensure_field(p,'noisefile','EsMatrixnoise_ltass.wav');    
        calibrationlevel=60;
        siglevel=50;
        profile='EsMatrix_ltass'; 
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri p.noisefile]; 
        warning('Fix noise as VlMatrix and EsMatrix')
    otherwise
        warning('Default speech material...')
        calibrationlevel=90; % play calibration signal at 90dBSPL
        siglevel=60;
        profile=sprintf('score_%d', p.calibration);
        cal_stim_id = 'gain0';
        filename=[p.speechmaterial_uri 'sin_1kHz_' num2str(siglevel-p.calibration) 'dBFS.wav'];
        warning('Fix noise as VlMatrix and EsMatrix')
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