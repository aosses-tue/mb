function Create_APEX_ITD_automated_learning
% function Create_APEX_ITD_automated_learning
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

% 'D:\Documenten-TUe\09-Training+activities\BEP-2015-2016\02-Eric\Meeting-20160316\APEX-experiments-handed-out\MATLAB\'
dir = 'D:\Documenten-TUe\09-Training+activities\BEP-2015-2016\02-Eric\20160322-from-Eric\';
filetemplate = [dir 'TEMPLATE.xml'];

p = []; % Empty struct

% Part 1: explaining how readfile_replace works:

p.thing_to_replace = 'this is what I want';
p.other_thing = 'this as well';

% Part 2: starting to automate:

p.procedures = []; % empty struct
p.datablock  = [];

dirname4stimuli = 'stimuli'; % Make sure inside this folder you put all stimuli
deviceid = 'wavdevice'; 

% L012-013  Each procedure
% A. Procedures
params = [];
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

type = 1; % adaptive
% procedure_id = {  'ProcedureLevel20', ... % id{1}
%                   'ProcedureLevel35'};    % id{2}

procedures_suffix = [20 35]; % test levels dB
procedures_id = {   sprintf('procedureLevel%.0f',procedures_suffix(1)), ...
                    sprintf('procedureLevel%.0f',procedures_suffix(2)) };

% % % One column for each procedure and amount of rows are the number of stimuli:
% %           For Procedure20 , Procedure35
% filename = {'wavefile-0.wav','wavefile-100.wav'; ...  % filename{1,1} = 'anyname';
%             'wavefile-1.wav','wavefile-101.wav'; ...
%             'wavefile-2.wav','wavefile-102.wav'; ...
%             'wavefile-3.wav','wavefile-103.wav'; ...
%             'wavefile-4.wav','wavefile-104.wav'; ...
%             'wavefile-5.wav','wavefile-105.wav'; ...
%             'wavefile-6.wav','wavefile-106.wav'; ...
%             'wavefile-7.wav','wavefile-107.wav'; ...
%             'wavefile-8.wav','wavefile-108.wav'; ...
%             'wavefile-9.wav','wavefile-109.wav'; ...
%             'wavefile-10.wav','wavefile-110.wav'; ...
%             'wavefile-11.wav','wavefile-111.wav'; ...
%             'wavefile-12.wav','wavefile-112.wav'; ...
%             'wavefile-13.wav','wavefile-113.wav'; ...
%             'wavefile-14.wav','wavefile-114.wav'; ...
%             'wavefile-15.wav','wavefile-115.wav'; ...
%             'wavefile-16.wav','wavefile-116.wav'};

for i = 1:17
    for j = 1:2
        filename{i,j} = sprintf('wavefile-%.0f.wav',i-1 + 100*(j-1) );
    end
end
nStimuli = size(filename,1); % each procedure

    
%     trials



%   A.1 Trials

% p.datablock = [];
 
% 
nTrials  =  2; % 2 ILD conditions
nStimuli_begin = [0 100]; % starts at   0 for procedure 20
                          % starts at 100 for procedure 35

for k = 1:length(procedures_id)
    
    p_il.trials = sprintf('<trials>\n');
    for i=1:nTrials
        trialid{i} = sprintf('trial%.0f-%.0f',procedures_suffix(k),i-1);
        answer  = 'button1'; % same 'answer' for every trial
        screen  = 'screen1'; % same 'screen'
        for j = 1:nStimuli
            stimulusid{j} = sprintf('stimulus%.0f-%.0f-%.0f', procedures_suffix(k), i-1, j-1);
            datablockid{j} = sprintf('datablock%.0f-%.0f-%.0f', procedures_suffix(k), i-1,  j-1); % same format as with stimulus (to be used in a3datablock)
            p.datablock       = [p.datablock a3datablock(datablockid{j},filename{j,k},deviceid)];
        end
        
        p_il.trials       = [p_il.trials a3trial(trialid{i}, screen, stimulusid, answer)];
%         
    end
    p_il.trials       = [p_il.trials sprintf('</trials>\n')];

    p.procedures = [sprintf('\n') p.procedures a3procedure(params, p_il.trials, type, procedures_id{k}) sprintf('\n')];
end

p.dirstimuli = dirname4stimuli;

XML_to_write = readfile_replace(filetemplate, p);

outputfile = [dir 'TEMPLATE-RESULT.xml']; % XML where I am going to write my output

fid     = fopen(outputfile, 'w');
fwrite(fid, XML_to_write);
fclose(fid);

disp('')
% XMLres = readfile_replace(filetemplate,p);
% 
% % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
