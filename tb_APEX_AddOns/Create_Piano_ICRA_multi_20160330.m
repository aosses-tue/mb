function Create_Piano_ICRA_multi_20160330(experiment_nr,do_skip)
% function Create_Piano_ICRA_multi_20160330(experiment_nr,do_skip)
%
% 1. Description:
%       Generate APEX experiments if do_skip is set to 0 (default) the audio
%       signals will be generated, otherwise only the APEX (*.apx) experiments.
%       The APEX experiments are based on the template 'piano_1_2_TEMPLATE.xml'
% 
% 2. Stand-alone example:
%   Experiments = [1 1.1 2 2.1  11 11.1 12 12.1 21 21.1 22 22.1];
%   Create_Piano_ICRA_multi_20160330(Experiments);
%   
%   Experiments = [1 1.1 2 2.1]; % only C2
%   do_skip = 0;
%   Create_Piano_ICRA_multi_20160330(Experiments,do_skip);
% 
%   Experiments = [11 11.1 12 12.1]; % only A4
%   do_skip = 0;
%   Create_Piano_ICRA_multi_20160330(Experiments,do_skip);
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also: r20160330_update_APEX_ICRA.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Adapted from: r20160330_update_APEX_ICRA.m
% Created on    : 25/03/2016
% Last update on: 25/03/2016 
% Last use on   : 25/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    do_skip = 1; % 1 = keeps already created audio's
end

dir_main   = [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments' delim 'Piano' delim 'Pilot-ICRA-v2' delim];
template_main           = [dir_main 'piano_multi_TEMPLATE.xml'];
template_1_proc         = [dir_main 'piano_multi_TEMPLATE_procedure.xml'];
template_2_datablocks   = [dir_main 'piano_multi_TEMPLATE_datablocks.xml'];
template_3_filters      = [dir_main 'piano_multi_TEMPLATE_filters.xml'];
template_4_stimuli      = [dir_main 'piano_multi_TEMPLATE_stimuli.xml'];
template_5_connections  = [dir_main 'piano_multi_TEMPLATE_connections.xml'];
template_6_cal          = [dir_main 'calibration-outs-only-L.xml'];

dir_out_XML     =  dir_main;

opts    = [];
p       = [];
p       = ef(p,'procedure',[]);
psub    = [];

%           0        1       2          3            4         5       6      7
pianos = {'GH05','GRAF28','JBS36','JBS51-4486','JBS51-4544','JBS50','JBS73','NS19'};
 
if nargin == 0
    experiment_nr = [1 1.1];
end

Nprocedures = length(experiment_nr);

fi = [];
fi_datablocks = [];
fi_pedestal = [];
rampout_len = 100; % ms
SNR4pede = 10; 

for nproc = 1:Nprocedures
    
    exp_nr = experiment_nr(nproc); % experiment being processed
    
    if experiment_nr < 10
        note_test = 'C2';
    elseif experiment_nr < 20
        note_test = 'A4';
    elseif experiment_nr < 30
        note_test = 'Csh5';
    end
     
    dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '05-loudness-balanced' delim note_test delim];
     
    switch exp_nr
        case 1 % C2 
            ID_1   = 4; % 4 = 'JBS51-4544'
            idx_1  = ID_1 + 1;
            take1  = [3 4]; 

            ID_2   = 7; % 7 = 'NS19';
            idx_2  = ID_2+1;
            take2     = [1 2]; 

        case 1.1
            ID_1   = 7; % 7 = 'NS19';
            idx_1  = ID_1+1;
            take1  = [1 2]; 

            ID_2   = 4; % 4 = 'JBS51-4544'
            idx_2  = ID_2 + 1;
            take2  = [3 4]; 

        case 2
            ID_1   = 1; % 1 = 'GRAF28'
            idx_1  = ID_1 + 1;
            take1  = [2 3]; 

            ID_2   = 7; % 7 = 'NS19';
            idx_2  = ID_2+1;
            take2     = [1 2]; 

        case 2.1
            ID_1   = 7; % 7 = 'NS19';
            idx_1  = ID_1+1;
            take1  = [1 2]; 

            ID_2   = 1; % 1 = 'GRAF28'
            idx_2  = ID_2 + 1;
            take2  = [2 3]; 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 11 % A4
            ID_1   = 4; % 4 = 'JBS51-4544'
            idx_1  = ID_1 + 1; 
            take1  = [3 4]; 

            ID_2   = 7; % 7 = 'NS19';
            idx_2  = ID_2+1; 
            take2     = [2 2]; 

        case 11.1
            ID_1   = 7; % 7 = 'NS19';
            idx_1  = ID_1+1; 
            take1  = [2 2]; 

            ID_2   = 4;
            idx_2  = ID_2 + 1; % 4 = 'JBS51-4544'
            take2  = [3 4]; 

        case 12
            ID_1   = 1; % 1 = 'GRAF28' (its take 1 was a pp)
            idx_1  = ID_1 + 1;
            take1  = [2 3]; 

            ID_2   = 7; % 7 = 'NS19';
            idx_2  = ID_2+1;
            take2     = [2 2]; 

        case 12.1
            ID_1   = 7; % 7 = 'NS19';
            idx_1  = ID_1+1;
            take1  = [2 2]; 

            ID_2   = 1;
            idx_2  = ID_2 + 1;
            take2  = [2 3]; 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 21 % Csh5
            ID_1   = 4;
            idx_1  = ID_1 + 1;
            take1  = [3 5]; 

            ID_2   = 7;
            idx_2  = ID_2+1;
            take2     = [2 4]; 

        case 21.1
            ID_1   = 7;
            idx_1  = ID_1+1;
            take1  = [2 4]; 

            ID_2   = 4;
            idx_2  = ID_2 + 1;
            take2  = [3 5]; 

        case 22
            ID_1   = 1;
            idx_1  = ID_1 + 1;
            take1  = [2 3]; 

            ID_2   = 7;
            idx_2  = ID_2+1;
            take2     = [2 4]; 

        case 22.1
            ID_1   = 7;
            idx_1  = ID_1+1;
            take1  = [2 4]; 

            ID_2   = 1;
            idx_2  = ID_2 + 1;
            take2  = [2 3]; 

    end
    
    piano1      = pianos{idx_1}; 
    piano2      = pianos{idx_2}; 
            
    dirstimuli = sprintf('Stimuli-%s',note_test);
    dir_out_stimuli = [dir_main dirstimuli delim];
    
    if do_skip == 0
        Mkdir(dir_out_XML);
        dir_out_aux = [dir_out_stimuli 'rawnoises' delim];
        Mkdir(dir_out_XML);
        Mkdir(dir_out_aux);
    end
     
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
    
        % 2. Aligning the sounds (looking at the maximum RMS value):
        [signal1, signal2, fs] = il_wavread(fname1, fname2, opts);
        noise1 = icra_noise4piano(signal1,fs);
        noise2 = icra_noise4piano(signal2,fs);
    
        noise3 = 0.5*noise1 + 0.5*noise2;
    
        out_tmp = auditoryfilterbank(noise3,fs);
        RMSrel  = rmsdb(out_tmp)-max(rmsdb(out_tmp));
        pede    = icra_noise4piano_pedestal(noise3,fs,RMSrel,0); % SNR of 0 dB, target SNR adjusted inside the experiment file
     
        % % 2.1. Saving ICRA noise:
        fname_out3p{i} = sprintf('pedestal-%s_%s-SNR-%.0f-dB.wav',fname1suffix,fname2suffix,0);

        signals{1+2*(i-1)} = signal1;
        signals{2+2*(i-1)} = signal2;
        noises{1+2*(i-1)}  = noise1;
        noises{2+2*(i-1)}  = noise2;

        fnames{1+2*(i-1)} = fname1suffix;
        fnames{2+2*(i-1)} = fname2suffix;

        if do_skip == 0
            Wavwrite(pede,fs,[dir_out_stimuli fname_out3p{i}]);
        end
         
    end
    
    L = length(signals{1}); % Length of the shortest audio segment
    for i = 1:length(signals)
        L = min(L,length(signals{i}));
    end

    if do_skip == 0
        for i = 1:length(signals)

            Wavwrite(signals{i}(1:L),fs,[dir_out_stimuli          fnames{i} '.wav']);
            Wavwrite( noises{i}(1:L),fs,[dir_out_aux     'noise-' fnames{i} '.wav']);

        end
    end

    for i = [1 3]
        for j = [2 4]
            noise3 = 0.5*noises{i}(1:L) + 0.5*noises{j}(1:L);
            if do_skip == 0
                Wavwrite(noise3,fs,sprintf('%snoise-%s_%s.wav',dir_out_stimuli,fnames{i},fnames{j}))
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
    
    ProcedureID{nproc}  = sprintf('%.0f%.0f',ID_1,ID_2);
    Piano1{nproc}       = sprintf('P%.0f',ID_1);
    Piano2{nproc}       = sprintf('P%.0f',ID_2);
    
    psub.ProcedureID    = ProcedureID{nproc};
    psub.test_note      = note_test;
    psub.Piano1label    = piano1;
    psub.Piano2label    = piano2;
    proc_tmp            = readfile_replace(template_1_proc,psub);
    
    p.procedure         = [p.procedure lf proc_tmp];
    
    psub = [];
    
    fi{end+1} = [fnames_Piano1{1} '.wav'];
    fi{end+1} = [fnames_Piano1{2} '.wav'];
    fi{end+1} = [fnames_Piano2{1} '.wav'];
    fi{end+1} = [fnames_Piano2{2} '.wav'];
    
    fi_stim{1,nproc} = fi{end-3};
    fi_stim{2,nproc} = fi{end-2}; 
    fi_stim{3,nproc} = fi{end-1};
    fi_stim{4,nproc} = fi{end};
    
    fi{end+1} = sprintf('noise-%s_%s.wav',fnames_Piano1{1},fnames_Piano2{1});
    fi{end+1} = sprintf('noise-%s_%s.wav',fnames_Piano1{1},fnames_Piano2{2});
    fi{end+1} = sprintf('noise-%s_%s.wav',fnames_Piano1{2},fnames_Piano2{1});
    fi{end+1} = sprintf('noise-%s_%s.wav',fnames_Piano1{2},fnames_Piano2{2});
    
    fi_noise{1,nproc}  = fi{end-3};
    fi_noise{2,nproc}  = fi{end-2}; % filenameNoise1_2
    fi_noise{3,nproc}  = fi{end-1}; % filenameNoise2_1
    fi_noise{4,nproc}  = fi{end}; % filenameNoise2_2
    
    fi{end+1}           = fname_out3p{1};
    fi_pedestal{end+1} = fi{end}; % p.filenamePedestal  = fi{end};

end

psub.datablock      = lf;
psub.dirstimuli     = dirstimuli;
psub.rampout_len    = rampout_len; % in ms

fi = unique(fi); % eliminating doubled file names
fi_pedestal = unique(fi_pedestal);

for i = 1:length(fi)
    
    suffixid = strsplit(fi{i},'.');
    dataid{i} = ['data_' suffixid{1}];
    
    psub.datablock = [psub.datablock a3datablock(dataid{i},fi{i},'wavdevice') lf];
    
end

p.datablocks      = readfile_replace(template_2_datablocks,psub);

%%%
psub = [];
psub.filters_pedestal = [];
psub.rampout_len = rampout_len;

for i = 1:length(fi_pedestal)
    
    suffixid    = strsplit(fi_pedestal{i},'.');
    dataid      = ['data_' suffixid{1}];
    filterid    = ['filter_' suffixid{1}];
    
    gain        = -SNR4pede;
    continuous  = 'true';
    psub.filters_pedestal = [psub.filters_pedestal lf a3dataloop(dataid,0,filterid,continuous,'wavdevice',gain)];
    
end

p.filters       = readfile_replace(template_3_filters,psub);

p.stimuli = [];

for nproc = 1:Nprocedures
    
    psub = [];
    
    psub.ProcedureID = ProcedureID{nproc};
    
    suffixid = strsplit(fi_stim{1,nproc},'.');
    psub.audio11 = ['data_' suffixid{1}];
    
    suffixid = strsplit(fi_noise{1,nproc},'.');
    psub.noise11 = ['data_' suffixid{1}];
    
    
    suffixid = strsplit(fi_stim{2,nproc},'.');
    psub.audio12 = ['data_' suffixid{1}];
    
    suffixid = strsplit(fi_noise{2,nproc},'.');
    psub.noise12 = ['data_' suffixid{1}];
    
    
    suffixid = strsplit(fi_stim{3,nproc},'.');
    psub.audio21 = ['data_' suffixid{1}];
    
    suffixid = strsplit(fi_noise{3,nproc},'.');
    psub.noise21 = ['data_' suffixid{1}];
    
    
    suffixid = strsplit(fi_stim{4,nproc},'.');
    psub.audio22 = ['data_' suffixid{1}];
    
    suffixid = strsplit(fi_noise{4,nproc},'.');
    psub.noise22 = ['data_' suffixid{1}];
    
    
    psub.Piano1 = Piano1{nproc};
    psub.Piano2 = Piano2{nproc};
    
    p.stimuli   = [p.stimuli lf readfile_replace(template_4_stimuli,psub)];
    
end

stimid = 'stimcaltone';
datablocks = psub.noise11; % arbitrary noise
calstim = a3stimulus(stimid, datablocks); %, fixedparameters, variableparameters, simultaneous);
    
p.stimuli = [p.stimuli lf calstim];
%%%

psub = [];
psub.connection = [];

fi_stim_unique = unique(fi_stim);

for i = 1:length(fi_stim_unique)
    
    stimid = strsplit(fi_stim_unique{i},'.');
    stimid = ['data_' stimid{1}];
    psub.connection = [psub.connection a3connection(stimid,0,'amp_use_test1',0)];
    
end

fi_noise_unique = unique(fi_noise);

for i = 1:length(fi_noise_unique)
    
    stimid = strsplit(fi_noise_unique{i},'.');
    stimid = ['data_' stimid{1}];
    psub.connection = [psub.connection a3connection(stimid,0,'rampout_inst_noise',0)];
    
end

fi_noise_unique = unique(fi_pedestal);

for i = 1:length(fi_noise_unique)
    
    stimid = strsplit(fi_noise_unique{i},'.');
    stimid = ['filter_' stimid{1}];
    psub.connection = [psub.connection a3connection(stimid,0,'amp_SNR',0)];
    
end

p.connections = readfile_replace(template_5_connections,psub);

p.calibration = readfile(template_6_cal);

%%%

XML_to_write = readfile_replace(template_main,p);
 
outputfile = [dir_main 'piano_multi_result-RENAME.xml']; % XML where I am going to write my output
fid     = fopen(outputfile, 'w');
fwrite(fid, XML_to_write);
fclose(fid);

disp('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    % then idx_1 is after
   [signal1 t1] = Do_alignment(t1,signal1,t1(samples_diff));
end

bAdditionalAlignment = 1;

if bAdditionalAlignment
    figure(100); 
    subplot(2,1,1)
    plot(t1,signal1, t2,signal2);
    haa = gca;
    
    [RMS_se1 t_se1] = rmsdb_sec(signal1,fs,window_length4max,0);
    [RMS_se2 t_se2] = rmsdb_sec(signal2,fs,window_length4max,0);
    
    subplot(2,1,2)
    plot(t_se1, RMS_se1, t_se2, RMS_se2);
    haa(end+1) = gca;
    linkaxes(haa,'x')
    legend('sig1','sig2')
    t1_diff = input('Enter the time t at which you want to synchronise for signal1 [s]: ' );
    t2_diff = input('Enter the time t at which you want to synchronise for signal2 [s]: ' );
    time_diff = t1_diff - t2_diff;    
    
    if time_diff ~= 0
        samples_diff = round(time_diff*fs);
        if time_diff < 0 
            % then idx_2 is after
           [signal2 t2] = Do_alignment(t2,signal2,t2(abs(samples_diff)));
        elseif time_diff > 0
            % then idx_1 is after
           [signal1 t1] = Do_alignment(t1,signal1,t1(samples_diff));
        end
        
    end
    
    figure(100); 
    subplot(2,1,1)
    plot(t1,signal1, t2,signal2);
    title('After manual alignment')
    
    [RMS_se1 t_se1] = rmsdb_sec(signal1,fs,window_length4max,0);
    [RMS_se2 t_se2] = rmsdb_sec(signal2,fs,window_length4max,0);
    
    subplot(2,1,2)
    plot(t_se1, RMS_se1, t_se2, RMS_se2);
    
end

[L idxL] = min([length(signal1) length(signal2)]);

RMS_se1 = rmsdb_sec(signal1,fs,window_length4cal);
RMS_se2 = rmsdb_sec(signal2,fs,window_length4cal);

max_1 = max(RMS_se1);
max_2 = max(RMS_se2);

delta_dB = max_2 - max_1;

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
