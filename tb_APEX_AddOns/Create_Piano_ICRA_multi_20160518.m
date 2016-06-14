function Create_Piano_ICRA_multi_20160518(pair_nr,do_skip,opts)
% function Create_Piano_ICRA_multi_20160518(experiment_nr,do_skip,opts)
%
% 1. Description:
%       Generate APEX experiments if do_skip is set to 0 (default) the audio
%       signals will be generated, otherwise only the APEX (*.apx) experiments.
%       The APEX experiments are based on the template 'piano_1_2_TEMPLATE.xml'
% 
% 2. Stand-alone example:
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

rampout_len = 0; % ms, stimuli already with ramp

bAdaptive = 1;
do_skip = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the experiments:
if bAdaptive
    nUp     = 1;
    nDown   = 2;
    max_value = 50;
    min_value = -25;
    start_value = 16; 
    stop_after  = 12;
    step1 = 4;
    step2 = 2; % 2
    step3 = 1; % 1
elseif bAdaptive == 0
    nPresentations = 16;
end

pause_between_stimuli = 200; % ms
intertrialscreen = 400; % ms

% Parameters of the stimuli:
i_noises      = 12; 
SNR4pede      = 40; 
gain4pede_ON  = 0;
gain4pede_OFF = -99;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    pair_nr = [1 2; 3 4];
end

if nargin < 2
    do_skip = 1; % 1 = keeps already created audio's
end

if nargin < 3
    opts = [];
end

opts = Ensure_field(opts,'dir_where',[Get_TUe_data_paths('piano') '04-PAPA' delim '06-Exp-TUe-1-similarity\02-to-be-used-in-AB-comparison' delim]);
opts = Ensure_field(opts,'dirstimuli','Stimuli-main-20160518');
opts = Ensure_field(opts,'SubjectID','SXX');
% opts = Ensure_field(opts,'RENAME',input('Input label for experiments (including quotes): '));

dir_where  = opts.dir_where;
RENAME     = opts.RENAME;
SubjectID  = opts.SubjectID;

dirstimuli = opts.dirstimuli; % name of the directory where the Stimuli will be stored

dir_main   = [Get_TUe_data_paths('ex_Experiments') '2015-APEX-my-experiments' delim 'Piano' delim 'Main-ICRA-v3' delim];

% folder_templates = 'Templates-v2';
folder_templates = 'Templates-v3'; % only minor changes with respect to v2

if bAdaptive == 1
    template_1_proc         = [dir_main folder_templates delim 'piano_multi_TEMPLATE_procedure.xml'];
    template_main           = [dir_main folder_templates delim 'piano_multi_TEMPLATE.xml'];
elseif bAdaptive == 0
    template_1_proc         = [dir_main folder_templates delim 'piano_multi_TEMPLATE_procedure-const.xml'];
    template_main           = [dir_main folder_templates delim 'piano_multi_TEMPLATE-interactive.xml'];
end

template_2_datablocks   = [dir_main folder_templates delim 'piano_multi_TEMPLATE_datablocks_one_lvl_up.xml'];
template_3_filters      = [dir_main folder_templates delim 'piano_multi_TEMPLATE_filters.xml'];
template_4_stimuli      = [dir_main folder_templates delim 'piano_multi_TEMPLATE_stimuli.xml'];
template_5_connections  = [dir_main folder_templates delim 'piano_multi_TEMPLATE_connections.xml'];
template_6_cal          = [dir_main folder_templates delim 'calibration-outs-only-L.xml'];

dir_out_XML     =  dir_main;
dir_out_stimuli = [dir_main dirstimuli delim];
Mkdir(dir_out_stimuli);

opts    = [];
p       = [];
p       = ef(p,'procedure',[]);
psub    = [];

%           0        1       2          3            4      5          6   
% pianos = {'GH05','GRAF28','JBS36','JBS51-4486','JBS51-4544','JBS73','NS19'};
files = Get_filenames(dir_where,'*.wav');

if do_skip == 0
    for i = 1:length(files)
        copyfile([dir_where files{i}],[dir_out_stimuli files{i}]);
        disp('')
    end
end

Nprocedures = size(pair_nr,1);

fi = [];
fi_datablocks = [];
fi_pedestal_nosort = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(['\tOutput dir for stimuli: %s\n' ...
         '\t %.0f noises will be generated \n' ...
         '\t Adaptive procedure\n'],dirstimuli,i_noises);
disp('Press any button to continue ')

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

for nproc = 1:Nprocedures
    
    ID = pair_nr(nproc,:); % experiment being processed
    idx     = ID;
    
    piano1 = strsplit(files{ID(1)},'-');
    if strcmp(piano1{1},'JBS51')
        piano1{1} = [piano1{1} piano1{2}];
    end
    piano1 = piano1{1}; 
    
    piano2 = strsplit(files{ID(2)},'-');
    if strcmp(piano2{1},'JBS51')
        piano2{1} = [piano2{1} piano2{2}];
    end
    piano2 = piano2{1}; 
    
    % 1. Sounds to be processed:

    fname1   = [dir_where files{ID(1)}];
    fname2   = [dir_where files{ID(2)}];

    if do_skip == 0
        [signal1, fs] = Wavread(fname1);
        [signal2, fs] = Wavread(fname2);

        for k = 1:i_noises
            noise1(:,k) = icra_noise4piano(signal1,fs);
            noise2(:,k) = icra_noise4piano(signal2,fs);
            noise_mix(:,k) = 0.5*noise1(:,k) + 0.5*noise2(:,k);
            disp(num2str(k))
        end

        out_tmp = auditoryfilterbank(noise_mix(:,1),fs);
        RMSrel  = rmsdb(out_tmp)-max(rmsdb(out_tmp));
        pede    = icra_noise4piano_pedestal(noise_mix(:,1),fs,RMSrel,0); % SNR of 0 dB, target SNR adjusted inside the experiment file

        signals{1} = signal1;
        signals{2} = signal2;
        noises{1}  = noise1(:,1);
        noises{2}  = noise2(:,1);
    end

     % 2.1. Saving pedestal noise:
    lvlPede = 60;
    fname_out3p{1} = sprintf('pedestal_%s_%s_%.0f_dB.wav',piano1,piano2,lvlPede);

    if do_skip == 0
        pede = setdbspl(pede,lvlPede);
        Wavwrite(pede,fs,[dir_out_stimuli fname_out3p{1}]);
    end
    
    % 2.2 Saving ICRA noise
    if do_skip == 0
        for k = 1:i_noises
            Wavwrite( noise_mix(:,k),fs,sprintf('%snoise_%s_%s_%.0f.wav',dir_out_stimuli,piano1,piano2,k));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fnames_Piano1 = [];
    fnames_Piano2 = [];
    idxtmp = 1; 
    for i = idxtmp;
        fnames_Piano1{end+1} = files{ID(1)};
    end
    idxtmp = 2; % 2 4
    for i = idxtmp;
        fnames_Piano2{end+1} = files{ID(2)};
    end
    
    ProcedureID{nproc}  = sprintf('%.0f%.0f',ID(1),ID(2));
    
    psub.ProcedureID    = ProcedureID{nproc};
    psub.Piano1label    = piano1;
    psub.Piano2label    = piano2;
    
    % bAdaptive procedure:
    if bAdaptive == 1
        psub.nUp            = nUp;
        psub.nDown          = nDown;
        psub.max_value      = max_value;
        psub.min_value      = min_value;
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
    
    fi{end+1} = [fnames_Piano1{1}];
    fi{end+1} = [fnames_Piano2{1}];
    
    fi_stim{1,nproc} = fi{end-1};
    fi_stim{2,nproc} = fi{end}; 
    
    fi{end+1} = sprintf('noise_%s_%s.wav',piano1,piano2);
     
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

    suffixid = strsplit(fi{i},'.'); % 4pedestal
    suffixid = strsplit(suffixid{1},'-');
    
    if strcmp(suffixid{1},'JBS51')
        suffixid{1} = [suffixid{1} suffixid{2}];
    end
    dataid{i} = ['data_' suffixid{1}];

    if ~strcmp(fi{i}(1:5),'noise')
        psub.datablock = [psub.datablock a3datablock(dataid{i},fi{i},'wavdevice') lf];
    else
        suffixid = suffixid{1};
        for k = 1:i_noises
            
            tmp_name = [suffixid '_' num2str(k)];
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
    suffixid    = strsplit(suffixid{1},'-');
    if strcmp(suffixid{1},'JBS51')
        suffixid{1} = [suffixid{1} suffixid{2}];
    end
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
    
    ID = pair_nr(nproc,:); % experiment being processed
    
    piano1 = strsplit(files{ID(1)},'-');
    if strcmp(piano1{1},'JBS51')
        piano1{1} = [piano1{1} piano1{2}];
    end
    piano1 = piano1{1}; 
    
    piano2 = strsplit(files{ID(2)},'-');
    if strcmp(piano2{1},'JBS51')
        piano2{1} = [piano2{1} piano2{2}];
    end
    piano2 = piano2{1}; 
    
    psub = [];
    psub.ProcedureID = ProcedureID{nproc};
    
    psub.audio1 = ['data_' piano1];
    psub.audio2 = ['data_' piano2];

    psub.noise_prefix = sprintf('noise_%s_%s',piano1,piano2);
    psub.pedestal_gains = '';
    
    for i = 1:length(fi_pedestal)
        suffixid    = strsplit(fi_pedestal{i},'.');
        if strcmp(filterid{i,2},psub.ProcedureID)
            psub.pedestal_gains = [psub.pedestal_gains '\n\t<parameter id="filter_' suffixid{1} '_gain">' num2str(gain4pede_ON) '</parameter>'];
        else
            psub.pedestal_gains = [psub.pedestal_gains '\n\t<parameter id="filter_' suffixid{1} '_gain">' num2str(gain4pede_OFF) '</parameter>'];
        end
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
    
    % stimid = strsplit(fi_stim_unique{i},'.');
    stimid = strsplit(fi_stim_unique{i},'-');
    if strcmp(stimid{1},'JBS51')
        stimid{1} = [stimid{1} stimid{2}];
    end
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
else
    p.interactiveSNR = '<entry type="double" description="SNR (dB)" expression="/apex:apex/filters[1]/filter[1]/gain[1]" default="20"/>';
end

p.intertrialscreen = intertrialscreen;
p.SubjectID = SubjectID;
XML_to_write = readfile_replace(template_main,p);
 
if length(RENAME) == 0
    RENAME = 'RENAME';
end

if bAdaptive == 1
    outputfile = [dir_main 'piano_multi-' RENAME '-adaptive' xmlextension]; % XML where I am going to write my output
elseif bAdaptive == 0
    outputfile = [dir_main 'piano_multi-' RENAME '-constant' xmlextension]; % XML where I am going to write my output
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
        outputfile = [dir_main 'piano_' RENAME '-const-' ProcedureID xmlextension]; % XML where I am going to write my output
    end
    
    fid     = fopen(outputfile, 'w');
    fwrite(fid, XML_to_write);
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
