function r20160329_update_APEX_ICRA
% function r20160329_update_APEX_ICRA
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 25/03/2016
% Last update on: 25/03/2016 
% Last use on   : 25/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_skip = 0; % 1 = keeps already created audio's

dir_main   = ['D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Hummer' delim];
dir_out = [dir_main 'ICRA' delim 'Stimuli' delim];
Mkdir(dir_out);

template = [dir_main 'Instrument-TEMPLATE-3AFC-running-m-AFC.xml.apx'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating the ICRA noise:
note_test = 'A4';
dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '03-Exported-as-segments' delim note_test delim];

opts = [];
    
piano_1 = 'JBS50';
take_1 = 3;

piano_2 = 'NS19';
take_2 = 1;
    
attenuate_by = 10; % Because these signals were measured very closed to the soundboard (they are very loud)
    
% 1. Sounds to be processed:
fname1suffix = [piano_1 '-' note_test '_' num2str(take_1)];
fname2suffix = [piano_2 '-' note_test '_' num2str(take_2)];
fname1   = [dir_where fname1suffix '.wav'];
fname2   = [dir_where fname2suffix '.wav'];
    
SNR4pede= 20;
if do_skip == 0
    % 2. Aligning the sounds (looking at the maximum RMS value):
    opts.attenuate_by = attenuate_by;
    [signal1, signal2, fs] = il_wavread(fname1, fname2, opts);
    noise1 = icra_noise4piano(signal1,fs);
    noise2 = icra_noise4piano(signal2,fs);

    noise3 = 0.5*noise1 + 0.5*noise2;

    out_tmp = auditoryfilterbank(noise3,fs);
    RMSrel  = rmsdb(out_tmp)-max(rmsdb(out_tmp));
    pede    = icra_noise4piano_pedestal(noise3,fs,RMSrel,SNR4pede);
end

% 2.1. Saving ICRA noise:
fname_out1 = sprintf('noise-%s.wav',fname1suffix);
fname_out2 = sprintf('noise-%s.wav',fname2suffix);
fname_out3 = sprintf('noise-%s-%s.wav',fname1suffix,fname2suffix);
fname_out3p = sprintf('pedestal-%s-%s-SNR-%.0f-dB.wav',fname1suffix,fname2suffix,SNR4pede);

if do_skip == 0
    Wavwrite(signal1,fs,[dir_out fname1suffix])
    Wavwrite(signal2,fs,[dir_out fname2suffix])

    Wavwrite(noise1,fs,[dir_out fname_out1]);
    Wavwrite(noise2,fs,[dir_out fname_out2]);
    Wavwrite(noise3,fs,[dir_out fname_out3]);
    Wavwrite(pede,fs,[dir_out fname_out3p]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I need to add roving

p = [];
p.wavepedestal      = fname_out3p;  % randomnoise-Fc-2510_BW-4980_SPL-77.wav
p.waveshapednoise   = fname_out3;   % hummer_noise_20151104.wav (make sure is 1 second long)
p.wavetest1         = [fname1suffix '.wav']; % meas-ac-4-dist-ane-HP+LP.wav
p.wavetest2         = [fname2suffix '.wav'];

XML_to_write = readfile_replace(template,p);

outputfile = [dir_out sprintf('piano-%s-%s.xml.apx',fname1suffix,fname2suffix)]; % XML where I am going to write my output
fid     = fopen(outputfile, 'w');
fwrite(fid, XML_to_write);
fclose(fid);
   
disp('')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal1, signal2, fs, sigorig1, sigorig2] = il_wavread(fname1,fname2,opts)
% function [signal1, signal2, fs, sigorig1, sigorig2] = il_wavread(fname1,fname2,opts)

if nargin < 3
    opts = [];
end
opts = Ensure_field(opts,'attenuate_by',0);
opts = Ensure_field(opts,'dur_ramp_down',150); % in ms, to be applied to truncated signal (the longer one)
opts = Ensure_field(opts,'window_length4max',10e-3); % in s
opts = Ensure_field(opts,'window_length4cal',0.5); % in s

attenuate_by = opts.attenuate_by;
dur_ramp_down = opts.dur_ramp_down;
window_length4max = opts.window_length4max;
window_length4cal = opts.window_length4cal;
% timeRMS2cal_bef = 150e-3; % assumes that the onset occurs at approx. t = 0.150 s
% timeRMS2cal_aft = 350e-3; % assumes that the target sound occurs in the next 0.350 s

[signal1, fs] = Wavread(fname1); 
[signal2, fs] = Wavread(fname2); 

t1 = ( 1:length(signal1) )/fs;
t2 = ( 1:length(signal2) )/fs;

[RMS_se1 t_se1] = rmsdb_sec(signal1,fs,window_length4max,0);
[RMS_se2 t_se2] = rmsdb_sec(signal2,fs,window_length4max,0);

%   2.1. Detecting the maximum and matching the levels:
[max_1 idx_1] = max(RMS_se1);
[max_2 idx_2] = max(RMS_se2);

idx_1 = find(t1 <= t_se1(idx_1),1,'last');
idx_2 = find(t2 <= t_se2(idx_2),1,'last');

samples_diff = abs(idx_2 - idx_1);
if idx_1 < idx_2 
    % then idx_2 is after
   [signal2 t2] = Do_alignment(t2,signal2,t2(samples_diff));
else
    % then idx_2 is after
   [signal1 t1] = Do_alignment(t1,signal1,t1(samples_diff));
end

[L idxL] = min([length(signal1) length(signal2)]);
% signal1 = Do_truncate(signal1,L);
% signal2 = Do_truncate(signal2,L);

RMS_se1 = rmsdb_sec(signal1,fs,window_length4cal);
RMS_se2 = rmsdb_sec(signal2,fs,window_length4cal);

max_1 = max(RMS_se1);
max_2 = max(RMS_se2);

delta_dB = max_2 - max_1;

if delta_dB >= 0 & delta_dB < 10 % then signal2 is louder (by less than 10 dB
    signal2 = From_dB(-abs(delta_dB)) * signal2; % attenuation of the louder signal
elseif abs(delta_dB) < 10
    signal1 = From_dB(-abs(delta_dB)) * signal1;
else
    error('The two signals you are trying to compare have very different levels')
end

signal1 = From_dB(-attenuate_by) * signal1;
signal2 = From_dB(-attenuate_by) * signal2;

if nargout > 3
    sigorig1 = signal1;
    sigorig2 = signal2;
end

signal1 = Do_truncate(signal1,L);
signal2 = Do_truncate(signal2,L);

if idxL == 1
    signal2 = Do_cos_ramp(signal2,fs,0,dur_ramp_down); % then signal2 was truncated
else
    signal1 = Do_cos_ramp(signal1,fs,0,dur_ramp_down); % then signal1 was truncated
end

