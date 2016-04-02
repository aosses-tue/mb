function Create_APEX_ITD_automated
% function Create_APEX_ITD_automated
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 23/03/2016
% Last update on: 23/03/2016 
% Last use on   : 23/03/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir = 'D:\Documenten-TUe\09-Training+activities\BEP-2015-2016\02-Eric\20160322-from-Eric\';
filetemplate = [dir 'TEMPLATE-selectrandom-ITD-variable-ILD-variable.xml-Ryan-Eric.apx'];
dirname4stimuli = 'stimuli'; % Make sure inside this folder you put all stimuli

deviceid = 'wavdevice';

% 1 column for each procedure
filename = {'wavefile-0.wav','wavefile-100.wav'; ...
            'wavefile-1.wav','wavefile-101.wav'; ...
            'wavefile-2.wav','wavefile-102.wav'; ...
            'wavefile-3.wav','wavefile-103.wav'; ...
            'wavefile-4.wav','wavefile-104.wav'; ...
            'wavefile-5.wav','wavefile-105.wav'; ...
            'wavefile-6.wav','wavefile-106.wav'; ...
            'wavefile-7.wav','wavefile-107.wav'; ...
            'wavefile-8.wav','wavefile-108.wav'; ...
            'wavefile-9.wav','wavefile-109.wav'; ...
            'wavefile-10.wav','wavefile-110.wav'; ...
            'wavefile-11.wav','wavefile-111.wav'; ...
            'wavefile-12.wav','wavefile-112.wav'; ...
            'wavefile-13.wav','wavefile-113.wav'; ...
            'wavefile-14.wav','wavefile-114.wav'; ...
            'wavefile-15.wav','wavefile-115.wav'; ...
            'wavefile-16.wav','wavefile-116.wav'};
nStimuli = size(filename,1); % each procedure

p = [];
p.procedures = [];
% % L001-011  Fixed

% % L012-013  Each procedure
% A. Procedures
params.presentations= 1;            % <presentations>1</presentations>
params.order        = 'sequential'; % <order>sequential</order>
params.choices      = 1;            %<choices>1</choices>
params.nUp          = 1;            % <nUp>1</nUp>
params.nDown        = 1;            % <nDown>1</nDown>
params.adapt_parameter = 'itd';     % <adapt_parameter>itd</adapt_parameter>
params.start_value  = 8;            % <start_value>8</start_value>
params.stop_after_type = 'reversals'; % <stop_after_type>reversals</stop_after_type>
params.stop_after   = 2;            % <stop_after>2</stop_after>
params.min_value    = -8;           % <min_value>-8</min_value>
params.max_value    = 8;            % <max_value>8</max_value>
params.larger_is_easier = 'false';  % <larger_is_easier>false</larger_is_easier>
params.repeat_first_until_correct = 'true'; % <repeat_first_until_correct>true</repeat_first_until_correct>
params.stepsizes    = sprintf('\n\t<stepsize begin="0" size="4"/>\n\t<stepsize begin="1" size="2"/>\n\t<stepsize begin="2" size="1"/>\n');

%   A.1 Trials
p_il.trials = sprintf('<trials>\n');
p.datablock = [];

procedures_suffix = [20 35];
procedures_id = {   sprintf('procedureLevel%.0f',procedures_suffix(1)), ...
                    sprintf('procedureLevel%.0f',procedures_suffix(2)) };

nTrials  =  1; % Tested up to now only with nTrial = 1
nStimuli_begin = [0 100]; % starts at   0 for procedure 20
                          % starts at 100 for procedure 35

for k = 1:length(procedures_id)
    for i=1:nTrials
        trialid{i} = sprintf('trial%.0f-%.0f',procedures_suffix(k),i-1);
        answer  = 'button1';
        screen  = 'screen1'; 
        for j = 1:nStimuli
            stimulusid{j} = sprintf('stimulus%.0f', nStimuli_begin(k)+j-1);
            datablockid{j} = sprintf('datablock%.0f', nStimuli_begin(k)+j-1); % same format as with stimulus (to be used in a3datablock)
            p.datablock       = [p.datablock a3datablock(datablockid{j},filename{j,k},deviceid)];
        end
        
        p_il.trials       = [p_il.trials a3trial(trialid{i}, screen, stimulusid, answer)];
        
    end
    p_il.trials       = [p_il.trials sprintf('</trials>\n')];

    p.procedures = [sprintf('\n') p.procedures a3procedure(params,p_il.trials,1,procedures_id{k}) sprintf('\n')];
end

p.dirstimuli = dirname4stimuli;

% B Datablocks:

% a3datablock(
% p.datablock = 

XMLres = readfile_replace(filetemplate,p);

outputfile = [dir 'my-first-experiment.xml'];
fid=fopen(outputfile, 'w');
fwrite(fid, XMLres);
fclose(fid);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
