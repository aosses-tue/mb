function r20160330_update_APEX_ICRA(experiment_nr,do_skip)
% function r20160330_update_APEX_ICRA(experiment_nr,do_skip)
%
% 1. Description:
%       Generate APEX experiments if do_skip is set to 0 (default) the audio
%       signals will be generated, otherwise only the APEX (*.apx) experiments.
%       The APEX experiments are based on the template 'piano_1_2_TEMPLATE.xml'
% 
% 2. Stand-alone example:
%   Experiments = [1 1.1 2 2.1  11 11.1 12 12.1 21 21.1 22 22.1];
%   for i = 1:length(Experiments)
%       r20160330_update_APEX_ICRA(Experiments(i));
%   end
%
%   Experiments = [1 1.1 2 2.1]; % only C2
%   do_skip = 0;
%   for i = 1:length(Experiments)
%       r20160330_update_APEX_ICRA(Experiments(i),do_skip);
%   end
% 
%   Experiments = [11 11.1 12 12.1]; % only A4
%   do_skip = 0;
%   for i = 1:length(Experiments)
%       r20160330_update_APEX_ICRA(Experiments(i),do_skip);
%   end
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 25/03/2016
% Last update on: 25/03/2016 
% Last use on   : 25/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    do_skip = 1; % 1 = keeps already created audio's
end

dir_main   = [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments' delim 'Piano' delim 'Pilot-ICRA-v2' delim];

% template = [dir_main 'Instrument-TEMPLATE-3AFC-running-m-AFC.xml.apx'];
template = [dir_main 'piano_1_2_TEMPLATE.xml'];

opts = [];
    
%           0        1       2          3            4         5       6      7
pianos = {'GH05','GRAF28','JBS36','JBS51-4486','JBS51-4544','JBS50','JBS73','NS19'};

% A4                2-3                             1-4                       2
% A4 (pp)            1                               -                        1

% Exp.  1-2 - A4
% Exp. 10   - Csh5

if nargin == 0
    experiment_nr = 12.1;
end

if experiment_nr < 10
    note_test = 'C2';
elseif experiment_nr < 20
    note_test = 'A4';
elseif experiment_nr < 30
    note_test = 'Csh5';
end

% dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '03-Exported-as-segments' delim note_test delim];
dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '05-loudness-balanced' delim note_test delim];

switch experiment_nr
    case 1 % C2
        ID_1   = 4;
        idx_1  = ID_1 + 1;
        piano1 = pianos{idx_1}; % 4 = 'JBS51-4544'
        take1  = [3 4]; 

        ID_2   = 7;
        idx_2  = ID_2+1;
        piano2 = pianos{idx_2}; % 7 = 'NS19';
        take2     = [1 2]; 

    case 1.1
        ID_1   = 7;
        idx_1  = ID_1+1;
        piano1 = pianos{idx_1}; % 7 = 'NS19';
        take1  = [1 2]; 

        ID_2   = 4;
        idx_2  = ID_2 + 1;
        piano2 = pianos{idx_2}; % 4 = 'JBS51-4544'
        take2  = [3 4]; 

    case 2
        ID_1   = 1;
        idx_1  = ID_1 + 1;
        piano1 = pianos{idx_1}; % 1 = 'GRAF28'
        take1  = [2 3]; 

        ID_2   = 7;
        idx_2  = ID_2+1;
        piano2 = pianos{idx_2}; % 7 = 'NS19';
        take2     = [1 2]; 

    case 2.1
        ID_1   = 7;
        idx_1  = ID_1+1;
        piano1 = pianos{idx_1}; % 7 = 'NS19';
        take1  = [1 2]; 

        ID_2   = 1;
        idx_2  = ID_2 + 1;
        piano2 = pianos{idx_2}; % 1 = 'GRAF28'
        take2  = [2 3]; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 11 % A4
        ID_1   = 4;
        idx_1  = ID_1 + 1;
        piano1 = pianos{idx_1}; % 4 = 'JBS51-4544'
        take1  = [3 4]; 

        ID_2   = 7;
        idx_2  = ID_2+1;
        piano2 = pianos{idx_2}; % 7 = 'NS19';
        take2     = [2 2]; 

    case 11.1
        ID_1   = 7;
        idx_1  = ID_1+1;
        piano1 = pianos{idx_1}; % 7 = 'NS19';
        take1  = [2 2]; 

        ID_2   = 4;
        idx_2  = ID_2 + 1;
        piano2 = pianos{idx_2}; % 4 = 'JBS51-4544'
        take2  = [3 4]; 

    case 12
        ID_1   = 1;
        idx_1  = ID_1 + 1;
        piano1 = pianos{idx_1}; % 1 = 'GRAF28' (its take 1 was a pp)
        take1  = [2 3]; 

        ID_2   = 7;
        idx_2  = ID_2+1;
        piano2 = pianos{idx_2}; % 7 = 'NS19';
        take2     = [2 2]; 

    case 12.1
        ID_1   = 7;
        idx_1  = ID_1+1;
        piano1 = pianos{idx_1}; % 7 = 'NS19';
        take1  = [2 2]; 

        ID_2   = 1;
        idx_2  = ID_2 + 1;
        piano2 = pianos{idx_2}; 
        take2  = [2 3]; 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 21 % Csh5
        ID_1   = 4;
        idx_1  = ID_1 + 1;
        piano1 = pianos{idx_1}; % 4 = JBS51-4544
        take1  = [3 5]; 

        ID_2   = 7;
        idx_2  = ID_2+1;
        piano2 = pianos{idx_2}; % 7 = 'NS19';
        take2     = [2 4]; 

    case 21.1
        ID_1   = 7;
        idx_1  = ID_1+1;
        piano1 = pianos{idx_1}; % 7 = 'NS19';
        take1  = [2 4]; 

        ID_2   = 4;
        idx_2  = ID_2 + 1;
        piano2 = pianos{idx_2}; % 4 = 'JBS51-4544'
        take2  = [3 5]; 

    case 22
        ID_1   = 1;
        idx_1  = ID_1 + 1;
        piano1 = pianos{idx_1}; % 1 = 'GRAF28'
        take1  = [2 3]; 

        ID_2   = 7;
        idx_2  = ID_2+1;
        piano2 = pianos{idx_2}; % 7 = 'NS19';
        take2     = [2 4]; 

    case 22.1
        ID_1   = 7;
        idx_1  = ID_1+1;
        piano1 = pianos{idx_1}; % 7 = 'NS19';
        take1  = [2 4]; 

        ID_2   = 1;
        idx_2  = ID_2 + 1;
        piano2 = pianos{idx_2}; 
        take2  = [2 3]; 

end

dirstimuli = sprintf('Stimuli-%s-P%.0f-P%.0f%s',note_test,ID_1,ID_2 );
dir_out = [dir_main 'Stage-3-loudness-balancing' delim dirstimuli delim];
dir_out_aux = [dir_out 'rawnoises' delim];
Mkdir(dir_out);
Mkdir(dir_out_aux);

fprintf('\tPiano 1 (test) = %s, takes %.0f, %.0f\n\tPiano 2 (ref ) = %s, takes %.0f, %.0f, is that right?\n', ...
            pianos{idx_1}, take1(1), take1(2), ...
            pianos{idx_2}, take2(1), take2(2));
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% pause()

for i = 1:length(take1)
    % 1. Sounds to be processed:
    fname1suffix = sprintf('P%.0ft%.0f',ID_1,take1(i)); % [ piano1 '-' note_test '_' num2str(take1)];
    fname2suffix = sprintf('P%.0ft%.0f',ID_2,take2(i));
    
    switch note_test
        case 'Csh5'
            note_test = 'Cd5';
    end
            
    fname1   = [dir_where piano1 '-' note_test '_' num2str(take1(i)) '.wav'];
    fname2   = [dir_where piano2 '-' note_test '_' num2str(take2(i)) '.wav'];

    SNR4pede= 10; 
    
    % 2. Aligning the sounds (looking at the maximum RMS value):
    [signal1, signal2, fs] = il_wavread(fname1, fname2, opts);
    noise1 = icra_noise4piano(signal1,fs);
    noise2 = icra_noise4piano(signal2,fs);

    noise3 = 0.5*noise1 + 0.5*noise2;

    out_tmp = auditoryfilterbank(noise3,fs);
    RMSrel  = rmsdb(out_tmp)-max(rmsdb(out_tmp));
    pede    = icra_noise4piano_pedestal(noise3,fs,RMSrel,SNR4pede);

    % % 2.1. Saving ICRA noise:
    fname_out3p{i} = sprintf('pedestal-%s_%s-SNR-%.0f-dB.wav',fname1suffix,fname2suffix,SNR4pede);

    signals{1+2*(i-1)} = signal1;
    signals{2+2*(i-1)} = signal2;
    noises{1+2*(i-1)}  = noise1;
    noises{2+2*(i-1)}  = noise2;
    
    fnames{1+2*(i-1)} = fname1suffix;
    fnames{2+2*(i-1)} = fname2suffix;
    
    if do_skip == 0
        Wavwrite(pede,fs,[dir_out fname_out3p{i}]);
    end
    
end

L = length(signals{1}); % Length of the shortest audio segment
for i = 1:length(signals)
    L = min(L,length(signals{i}));
end

if do_skip == 0
    for i = 1:length(signals)

        Wavwrite(signals{i}(1:L),fs,[dir_out              fnames{i} '.wav']);
        Wavwrite( noises{i}(1:L),fs,[dir_out_aux 'noise-' fnames{i} '.wav']);

    end
end

for i = [1 3]
    for j = [2 4]
        noise3 = 0.5*noises{i}(1:L) + 0.5*noises{j}(1:L);
        if do_skip == 0
            Wavwrite(noise3,fs,sprintf('%snoise-%s_%s.wav',dir_out,fnames{i},fnames{j}))
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames_Piano1 = [];
fnames_Piano2 = [];
idxtmp = [1 3];
for i = idxtmp;
    fnames_Piano1{end+1} = fnames{i};
end
idxtmp = [2 4];
for i = idxtmp;
    fnames_Piano2{end+1} = fnames{i};
end

p = [];
p.ProcedureID       = sprintf('%.0f%.0f',ID_1,ID_2);
p.Piano1            = sprintf('P%.0f',ID_1);
p.Piano2            = sprintf('P%.0f',ID_2);
p.rampout_len       = 100; % ms 

p.filenamePiano1_1  = [fnames_Piano1{1} '.wav'];
p.filenamePiano1_2  = [fnames_Piano1{2} '.wav']; 
p.filenamePiano2_1  = [fnames_Piano2{1} '.wav'];
p.filenamePiano2_2  = [fnames_Piano2{2} '.wav'];

p.filenameNoise1_1  = sprintf('noise-%s_%s.wav',fnames_Piano1{1},fnames_Piano2{1}); 
p.filenameNoise1_2  = sprintf('noise-%s_%s.wav',fnames_Piano1{1},fnames_Piano2{2});
p.filenameNoise2_1  = sprintf('noise-%s_%s.wav',fnames_Piano1{2},fnames_Piano2{1});
p.filenameNoise2_2  = sprintf('noise-%s_%s.wav',fnames_Piano1{2},fnames_Piano2{2});

p.filenamePedestal  = fname_out3p{1};  
p.dirstimuli        = dirstimuli;

XML_to_write = readfile_replace(template,p);

outputfile = [dir_out sprintf('piano-P%.0f-P%.0f-%s.xml.apx',ID_1,ID_2,note_test)]; % XML where I am going to write my output
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
opts = Ensure_field(opts,'attenuate_by',[0 0]);
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
elseif samples_diff ~= 0
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

% if delta_dB >= 0 & delta_dB < 20 % then signal2 is louder (by less than 10 dB
%     signal2 = From_dB(-abs(delta_dB)) * signal2; % attenuation of the louder signal
% elseif abs(delta_dB) < 20
%     signal1 = From_dB(-abs(delta_dB)) * signal1;
% else
%     error('The two signals you are trying to compare have very different levels')
% end

signal1 = From_dB(-attenuate_by(1)) * signal1;
signal2 = From_dB(-attenuate_by(2)) * signal2;

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
