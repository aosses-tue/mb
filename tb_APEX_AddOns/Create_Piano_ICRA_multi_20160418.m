function Create_Piano_ICRA_multi_20160418(experiment_nr,do_skip,bAdaptive,opts)
% function Create_Piano_ICRA_multi_20160418(experiment_nr,do_skip,bAdaptive,opts)
%
% 1. Description:
%       Generate APEX experiments if do_skip is set to 0 (default) the audio
%       signals will be generated, otherwise only the APEX (*.apx) experiments.
%       The APEX experiments are based on the template 'piano_1_2_TEMPLATE.xml'
% 
% 2. Stand-alone example:
%   Experiments = [1 1.1 2 2.1  11 11.1 12 12.1 21 21.1 22 22.1];
%   Create_Piano_ICRA_multi_20160418(Experiments);
%   
%   Experiments = [1 1.1 2 2.1]; % only C2
%   do_skip = 0;
%   Create_Piano_ICRA_multi_20160418(Experiments,do_skip);
% 
%   Experiments = 17:0.1:17.5;
%   do_skip = 0;
%   Create_Piano_ICRA_multi_20160418(Experiments,do_skip);
% 
%   Experiments = [17 17.1]; 
%   do_skip = 0;
%   Create_Piano_ICRA_multi_20160418(Experiments,do_skip);
% 
% 3. Additional info:
%       Tested cross-platform: No
%       See also: r20160330_update_APEX_ICRA, r20160406_update_APEX_ICRA
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Adapted from: r20160406_update_APEX_ICRA.m
% Created on    : 25/03/2016
% Last update on: 25/03/2016 
% Last use on   : 25/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    do_skip = 1; % 1 = keeps already created audio's
end

if nargin < 3
    bAdaptive = 1;
end

if nargin == 0
    experiment_nr = [17 17.1];
end

if nargin < 4
    opts = [];
end

opts = Ensure_field(opts,'dirstimuli','Stimuli-pilot-20160502');
dirstimuli = opts.dirstimuli; % name of the directory where the Stimuli will be stored

dir_main   = [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments' delim 'Piano' delim 'Pilot-ICRA-v2' delim];
template_main           = [dir_main 'Templates-v1' delim 'piano_multi_TEMPLATE.xml'];
if bAdaptive == 1
    template_1_proc         = [dir_main 'Templates-v1' delim 'piano_multi_TEMPLATE_procedure.xml'];
elseif bAdaptive == 0
    template_1_proc         = [dir_main 'Templates-v1' delim 'piano_multi_TEMPLATE_procedure-const.xml'];
end

template_2_datablocks   = [dir_main 'Templates-v1' delim 'piano_multi_TEMPLATE_datablocks.xml'];
template_3_filters      = [dir_main 'Templates-v1' delim 'piano_multi_TEMPLATE_filters.xml'];
template_4_stimuli      = [dir_main 'Templates-v1' delim 'piano_multi_TEMPLATE_stimuli.xml'];
template_5_connections  = [dir_main 'Templates-v1' delim 'piano_multi_TEMPLATE_connections.xml'];
template_6_cal          = [dir_main 'Templates-v1' delim 'calibration-outs-only-L.xml'];

dir_out_XML     =  dir_main;

opts    = [];
p       = [];
p       = ef(p,'procedure',[]);
psub    = [];

%           0        1       2          3            4         5       6      7
pianos = {'GH05','GRAF28','JBS36','JBS51-4486','JBS51-4544','JBS50','JBS73','NS19'};
 
Nprocedures = length(experiment_nr);

fi = [];
fi_datablocks = [];
fi_pedestal_nosort = [];
rampout_len = 100; % ms
SNR4pede    = 20; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the experiments:
nUp     = 1;
nDown   = 2;
max_value = 50;
start_value = 8; % 16
stop_after = 12;
step1 = 4;
step2 = 4; % 2
step3 = 2; % 1
pause_between_stimuli = 250; % ms
nPresentations = 12;

% Parameters of the stimuli:
bL = 1; % if 1 = then the manual set L-length is going to be used.
fs_theo    = 44100;
dur_piano_samples = 0.1+1+rampout_len/1000;

L          = round(dur_piano_samples*fs_theo); 
loud2use   = 18; % sones
suffixloud = sprintf('-%.0f-sone',loud2use);
i_noises   = 12;
gain4pede_ON = 0;
gain4pede_OFF = -99;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['\tStimuli will be %.0f [ms] long\n ' ...
         '\tOutput dir for stimuli: %s\n' ...
         '\t %.0f noises will be generated \n' ...
         '\t Experimental procedure = %.0f (1 = adaptive; 0 = constant)\n'],L/fs_theo*1000,dirstimuli,i_noises,bAdaptive);
disp('Press any button to continue ')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
pause()

for nproc = 1:Nprocedures
    
    exp_nr = experiment_nr(nproc); % experiment being processed
    
    if exp_nr < 10
        if exp_nr < 5
            note_test = 'C2';
        else
            note_test = 'F3';
        end
        
    elseif exp_nr < 20
        if exp_nr < 15
            note_test = 'A4';
        else
            note_test = 'C4';
        end
        
    elseif exp_nr < 30
        note_test = 'Csh5';
    end
     
    dir_where = [Get_TUe_data_paths('piano') '04-PAPA' delim '05-loudness-balanced' delim note_test delim];
     
    switch exp_nr
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 3 % C2
            ID   = [0 2]; % Piano 0, piano 2
            take = [1 1]; % take  1, take  1
            
        case 3.1
            ID   = [2 4]; 
            take = [1 3]; 
        
        case 3.2
            ID   = [4 0]; 
            take = [3 1]; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 16 % C4   
            ID   = [4 4]; % Piano 4, piano 4 (JBS51-4544)
            take = [2 3]; % take  2, take  3
            
        case 16.1  
            ID   = [4 4]; 
            take = [1 3];

        case 16.2 
            ID   = [4 4];
            take = [4 2];
            
        case 16.4 
            ID   = [4 4]; 
            take = [3 2]; % JBS51-4544
            
        case 17
            ID   = [2 2]; % JBS36, intra piano comparison
            take = [1 2];
            
        case 17.1
            ID   = [2 2]; % JBS36, intra piano comparison
            take = [2 1];
        
        case 17.2 
            ID   = [2 4]; % JBS36, JBS51-4544, inter piano, similar
            take = [2 3];
        
        case 17.3
            ID   = [4 2]; % JBS51-4544, JBS36, inter piano, similar
            take = [3 2];
            
        case 17.4
            ID   = [0 2]; % GH05, JBS36, inter piano, dissimilar
            take = [2 2];
        
        case 17.5
            ID   = [2 0]; % JBS36, GH05, inter piano, dissimilar
            take = [2 2];
            
        case 17.6
            ID   = [0 4]; % GH05, JBS51-4544
            take = [2 3];
            
        case 17.7
            ID   = [4 0]; % JBS51-4544, GH05
            take = [3 2];
            
        case 17.8
            ID   = [4 6]; 
            take = [4 3];
            
        case 18 % C4
            ID   = [0 2]; 
            take = [2 2];
            
        case 18.1
            ID   = [2 4]; 
            take = [2 3];

        case 18.2
            ID   = [4 0]; 
            take = [3 2];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 25 % Csh5
            ID   = [0 2]; 
            take = [1 3];
            
        case 25.1
            ID   = [2 4]; 
            take = [3 3];

        case 25.2
            ID   = [4 0]; 
            take = [3 1];

    end
    idx = ID + 1;
    
    piano1      = pianos{idx(1)}; 
    piano2      = pianos{idx(2)}; 
            
    dir_out_stimuli = [dir_main dirstimuli delim];
    
    if do_skip == 0
        Mkdir(dir_out_XML);
        dir_out_aux = [dir_out_stimuli 'rawnoises' delim];
        Mkdir(dir_out_XML);
        Mkdir(dir_out_aux);
    end
     
    for i = 1:length(take(1))
        % 1. Sounds to be processed:
        fname1suffix = sprintf('P%.0ft%.0f',ID(1),take(1)); % [ piano1 '-' note_test '_' num2str(take1)];
        fname2suffix = sprintf('P%.0ft%.0f',ID(2),take(2));
        
        switch note_test
            case 'Csh5'
                note_test = 'Cd5';
        end
                
        fname1   = [dir_where piano1 '-' note_test '_' num2str(take(1)) suffixloud '.wav'];
        fname2   = [dir_where piano2 '-' note_test '_' num2str(take(2)) suffixloud '.wav'];
    
        if do_skip == 0
            % 2. Aligning the sounds (looking at the maximum RMS value):
            [signal1, signal2, fs] = il_wavread(fname1, fname2, opts);
            
            if bL == 0
                L = min( length(signal1),length(signal2) ); % Length of the shortest audio segment
            end
            
            for k = 1:i_noises
                noise1(:,k) = icra_noise4piano(signal1(1:L,:),fs);
                noise2(:,k) = icra_noise4piano(signal2(1:L,:),fs);
                noise3(:,k) = 0.5*noise1(:,k) + 0.5*noise2(:,k);
                disp(num2str(k))
            end

            out_tmp = auditoryfilterbank(noise3(:,1),fs);
            RMSrel  = rmsdb(out_tmp)-max(rmsdb(out_tmp));
            pede    = icra_noise4piano_pedestal(noise3(:,1),fs,RMSrel,0); % SNR of 0 dB, target SNR adjusted inside the experiment file
     
            signals{1+2*(i-1)} = signal1;
            signals{2+2*(i-1)} = signal2;
            noises{1+2*(i-1)}  = noise1(:,1);
            noises{2+2*(i-1)}  = noise2(:,1);
        end
        
        fnames{1+2*(i-1)} = [note_test fname1suffix];
        fnames{2+2*(i-1)} = [note_test fname2suffix];
                
        % % 2.1. Saving ICRA noise:
        lvlPede = 60;
        
        % fname_out3p{i} = sprintf('%s_pedestal_%s_%s_SNR_%.0f_dB.wav',note_test,fname1suffix,fname2suffix,0);
        fname_out3p{i} = sprintf('%s_pedestal_%s_%s_%.0f_dB.wav',note_test,fname1suffix,fname2suffix,lvlPede);
        
        if do_skip == 0
            pede = setdbspl(pede,lvlPede);
            Wavwrite(pede,fs,[dir_out_stimuli fname_out3p{i}]);
        end
         
    end
    
    if do_skip == 0
        
        for i = 1:length(signals)

            Wavwrite(signals{i}(1:L),fs,[dir_out_stimuli          fnames{i} '.wav']);
            Wavwrite( noises{i}(1:L),fs,[dir_out_aux     'noise_' fnames{i} '.wav']);
            
        end
        
        for k = 1:i_noises
            Wavwrite( noise3(1:L,k),fs,sprintf('%snoise_%s_%s_%.0f.wav',dir_out_stimuli,fnames{1},fnames{2},k));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fnames_Piano1 = [];
    fnames_Piano2 = [];
    idxtmp = 1; % 1 3
    for i = idxtmp;
        fnames_Piano1{end+1} = fnames{i};
    end
    idxtmp = 2; % 2 4
    for i = idxtmp;
        fnames_Piano2{end+1} = fnames{i};
    end
    
    if ID(1) == ID(2)
        ProcedureID{nproc}  = sprintf('%.0f_t%.0f_t%.0f',ID(1),take(1),take(2));
    else
        ProcedureID{nproc}  = sprintf('%.0f%.0f',ID(1),ID(2));
    end
    Piano1{nproc}       = sprintf('P%.0f',ID(1));
    Piano2{nproc}       = sprintf('P%.0f',ID(2));
    
    psub.ProcedureID    = ProcedureID{nproc};
    psub.test_note      = note_test;
    psub.Piano1label    = piano1;
    psub.Piano2label    = piano2;
    if bAdaptive == 1
        psub.nUp            = nUp;
        psub.nDown          = nDown;
        psub.max_value      = max_value;
        psub.stop_after     = stop_after;
        psub.start_value    = start_value;
        psub.step1    = step1;
        psub.step2    = step2;
        psub.step3    = step3;
    elseif bAdaptive == 0
        psub.nPresentations = nPresentations;
    end
    psub.pause_between_stimuli = pause_between_stimuli;
    
    proc_tmp            = readfile_replace(template_1_proc,psub);
    
    p.procedure         = [p.procedure lf proc_tmp];
    
    exp2eval = sprintf('p_this%.0f.procedure = proc_tmp; % to generate individual experiments',nproc);
    eval(exp2eval);
    exp2eval = sprintf('p_this%.0f.ProcedureID = ProcedureID{nproc};',nproc);
    eval(exp2eval);
    
    psub = [];
    
    fi{end+1} = [fnames_Piano1{1} '.wav'];
    fi{end+1} = [fnames_Piano2{1} '.wav'];
    
    fi_stim{1,nproc} = fi{end-1};
    fi_stim{2,nproc} = fi{end}; 
    
    fi{end+1} = sprintf('noise_%s_%s.wav',fnames_Piano1{1},fnames_Piano2{1});
     
    % fi_noise{1,nproc}  = fi{end};
    
    fi{end+1}           = fname_out3p{1};
    fi_pedestal_nosort{end+1} = fi{end}; % p.filenamePedestal  = fi{end};

end

fi_noise = '';

psub.datablock      = lf;
psub.dirstimuli     = dirstimuli;
psub.rampout_len    = rampout_len; % in ms

fi = unique(fi); % eliminating doubled file names
fi_pedestal = unique(fi_pedestal_nosort);

for i = 1:length(fi)

    suffixid = strsplit(fi{i},'.');
    dataid{i} = ['data_' suffixid{1}];

    if ~strcmp(fi{i}(1:5),'noise')
        psub.datablock = [psub.datablock a3datablock(dataid{i},fi{i},'wavdevice') lf];
    else
        % dataid{i} = ['data_' suffixid{1}];
        for k = 1:i_noises
            tmp_name = [suffixid{1} '_' num2str(k)];
            fi{end+1} = [tmp_name '.wav'];
            if k == 1
                fi{i} = fi{end};
                fi_noise{end+1} = fi{i};
                fi(end) = []; 
                % fi_noise{end-1} = [];
            else
                fi_noise{end+1} = fi{end};
            end
            psub.datablock = [psub.datablock a3datablock(['data_' tmp_name],[tmp_name '.wav'],'wavdevice') lf];
        end
    end

end

p.datablocks      = readfile_replace(template_2_datablocks,psub);

%%%
psub = [];
psub.filters_pedestal = [];
psub.rampout_len = rampout_len;

for i = 1:length(fi_pedestal)
    
    suffixid    = strsplit(fi_pedestal{i},'.');
    dataid      = ['data_' suffixid{1}];
    filterid{i,1} = suffixid{1};
    filterid{i,2} = ProcedureID{i};
    gain        = 0;
    continuous  = 'true';
    psub.filters_pedestal = [psub.filters_pedestal lf a3dataloop(dataid,0,['filter_' filterid{i,1}],continuous,'wavdevice',gain)];
    
end

psub.SNR4pede   = SNR4pede;
p.filters       = readfile_replace(template_3_filters,psub);

p.stimuli = [];

for nproc = 1:Nprocedures
    
    psub = [];
    
    psub.ProcedureID = ProcedureID{nproc};
    
    suffixid = strsplit(fi_stim{1,nproc},'.');
    psub.audio1 = ['data_' suffixid{1}];
    
    suffixid = strsplit(fi_stim{2,nproc},'.');
    psub.audio2 = ['data_' suffixid{1}];

    psub.noise_prefix = sprintf('noise_%s_%s',fnames{1},fnames{2});
    psub.pedestal_gains = '';
    
    for i = 1:length(fi_pedestal)
        suffixid    = strsplit(fi_pedestal{i},'.');
        if strcmp(filterid{i,2},psub.ProcedureID)
            psub.pedestal_gains = [psub.pedestal_gains '\n\t<parameter id="filter_' suffixid{1} '_gain">' num2str(gain4pede_ON) '</parameter>'];
        else
            psub.pedestal_gains = [psub.pedestal_gains '\n\t<parameter id="filter_' suffixid{1} '_gain">' num2str(gain4pede_OFF) '</parameter>'];
        end
        warning('continuar aqui')
    end
    p.stimuli   = [p.stimuli lf readfile_replace(template_4_stimuli,psub)];
    
end

stimid = 'stimcaltone';
datablocks = ['data_' psub.noise_prefix '_1']; % arbitrary noise
calstim = a3stimulus(stimid, datablocks); %, fixedparameters, variableparameters, simultaneous);
    
p.stimuli = [p.stimuli lf calstim];
%%%

psub = [];
psub.connection = [];

fi_stim_unique = unique(fi_stim);

for i = 1:length(fi_stim_unique)
    
    stimid = strsplit(fi_stim_unique{i},'.');
    stimid = ['data_' stimid{1}];
    psub.connection = [psub.connection a3connection(stimid,0,'rampout',0)];
    
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
    psub.connection = [psub.connection a3connection(stimid,0,'pede_SNR',0)];
    
end

p.connections = readfile_replace(template_5_connections,psub);

p.calibration = readfile(template_6_cal);

xmlextension = '.xml.apx';
%%%
if bAdaptive == 1
    p.interactiveSNR = '';
elseif bAdaptive == 0
    p.interactiveSNR = '<entry type="double" description="SNR (dB)" expression="/apex:apex/filters[1]/filter[1]/gain[1]" default="20"/>';
end
XML_to_write = readfile_replace(template_main,p);
 
RENAME = input('Input label for experiments (including quotes): ');
if length(RENAME) == 0
    RENAME = 'RENAME';
end

if bAdaptive == 1
    outputfile = [dir_main 'piano_multi-' RENAME '-adaptive' xmlextension]; % XML where I am going to write my output
elseif bAdaptive == 0
    outputfile = [dir_main 'piano_multi-' RENAME '-constant' xmlextension]; 
end
    
fid     = fopen(outputfile, 'w');
fwrite(fid, XML_to_write);
fclose(fid);

for nproc = 1:Nprocedures
    p_this = p;
    exp2eval = sprintf('p_this.procedure = p_this%.0f.procedure;',nproc);
    eval(exp2eval);
    exp2eval = sprintf('ProcedureID = p_this%.0f.ProcedureID;',nproc);
    eval(exp2eval);
    
    XML_to_write = readfile_replace(template_main,p_this);
 
    if bAdaptive == 1
        outputfile = [dir_main 'piano_' RENAME '-adaptive-' ProcedureID xmlextension]; % XML where I am going to write my output
    elseif bAdaptive == 0
        outputfile = [dir_main 'piano_' RENAME '-constant-' ProcedureID xmlextension]; 
    end

    fid     = fopen(outputfile, 'w');
    fwrite(fid, XML_to_write);
    fclose(fid);
end

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

bAdditionalAlignment = 0;

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
    title('After manual alignment (zero padding will be added to both signals)')
    
    dur_ramp_up = 10; % 10-ms
    signal1 = Do_cos_ramp(signal1,fs,dur_ramp_up, 0); % in case of noisy recordings (to avoid clicks)
    signal2 = Do_cos_ramp(signal2,fs,dur_ramp_up, 0); % in case of noisy recordings 
    sil2add = Gen_silence(abs(time_diff),fs); 
    signal1 = [sil2add ; signal1];
    signal2 = [sil2add ; signal2];
        
    [RMS_se1 t_se1] = rmsdb_sec(signal1,fs,window_length4max,0);
    [RMS_se2 t_se2] = rmsdb_sec(signal2,fs,window_length4max,0);
    
    subplot(2,1,2)
    plot(t_se1, RMS_se1, t_se2, RMS_se2);
    
else
    warning('Additional alignment turned off, be aware of this')
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

signal2 = Do_cos_ramp(signal2,fs,0,dur_ramp_down); % then signal2 was truncated
signal1 = Do_cos_ramp(signal1,fs,0,dur_ramp_down); % then signal1 was truncated
