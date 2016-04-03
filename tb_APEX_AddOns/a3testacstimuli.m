function a3testacstimuli(savefilename, path, ncols, opts)
% function a3testacstimuli(savefilename, path, ncols, opts)
% 
% 1. Description:
%       Creates an APEX3 experiment file that can present any of the .wav
%       files under the directory 'paths'.
%       'ncols' correspond to the number of columns to be displayed. If there
%       are 5 wav files in path, then buttons are going to be displayed in a
%       2 x 3 matrix
%       It uses the XML template file 'a3testacstimuli.xml' and it requires
%       the APEX MATLAB toolbox
% 
% 2. Stand-alone example (Windows example):
%   % 2.1 Windows example (no calibration tone):
%       path = 'D:\Output\Daniel1997_test_20141126\'; 
%       savefilename = 'test-XML.xml';
%       a3testacstimuli(savefilename, path,3);
% 
%   % 2.2 Windows example (calibration tone):
%       path = 'D:\Output\Daniel1997_test_20141126\'; 
%       savefilename = 'test-XML.xml';
%       opts.bCalTone = 1;
%       a3testacstimuli(savefilename, path,3,opts);
% 
%   % 2.3 Windows example (calibration tone):
%       path = 'D:\Output\Daniel1997_test_20141126\'; 
%       savefilename = 'test-XML.xml';
%       opts.bCalTone = 1;
%       opts.presentation = 'diotic';
%       a3testacstimuli(savefilename, path,3,opts);
% 
% Programmed by the APEX3 team
% Created on    : 2013-2014
% Adapted by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Last update on: 11/06/2015 
% Last use on   : 02/04/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    opts = [];
end

if nargin < 3
    ncols = 3;
end

if nargin < 2
    try
        path = uigetdir(Get_TUe_paths('outputs'),'Please choose a folder containing one or more wav files');
    catch
        path = uigetdir(pwd,'Please choose a folder containing one or more wav files');
    end
    path = [path delim];
end

opts = Ensure_field(opts,'bCalTone',0);
opts = Ensure_field(opts,'presentation','monaural');
bCalTone = opts.bCalTone;
presentation = opts.presentation;

experiment.trials='';
experiment.screen='';
experiment.datablocks='';
experiment.map='';
experiment.stimuli='';

try 
    dir_templates = Get_TUe_paths('tb_APEX_AddOns_templates');
catch
    dir_templates = '';
end

try
    path=makedirend(path);
end

files=dir([path '*.wav']);
buttons=cell(length(files),1);

for f=1:length(files)
    filename=files(f).name;
    filename_noext=filename(1:end-4);
    
    datablock=['d' filename_noext];
    stimulus=['s' filename_noext];
    trial=['t' filename_noext];
    button=['button' filename_noext];
    experiment.datablocks=[experiment.datablocks a3datablock(datablock, filename, 'wavdevice')]; % 'soundcard' %for old templates
    experiment.trials=[experiment.trials a3trial(trial, 'screen', stimulus, button)];
    experiment.stimuli=[experiment.stimuli a3stimulus(stimulus, datablock)];
    
    buttons{f}=filename_noext;
end

if bCalTone == 0
    experiment.calstimulus = a3stimulus_empty('calstimulus');
else
    disp('Select a file to use as the calibration tone: ')
    for i = 1:length(files)
        fprintf('(Press %.0f) %s\n',i,files(i).name);
    end
    idxcal = input('Enter your choice: ');
    
    filename=files(idxcal).name;
    filename_noext=filename(1:end-4);
    datablock=['d' filename_noext];
    
    experiment.calstimulus = a3stimulus('calstimulus', datablock);
end

nrows = ceil(length(buttons)/ncols);

buttons=[buttons; cell(nrows*ncols-length(buttons),1)];

buttons = transpose( reshape(buttons, ncols, nrows) ); % trick to get [1 2 3; 4 5 6] instead of [1 3 5;2 4 6]

experiment.buttonlayout = a3buttonlayout(buttons, 1,2);
experiment.buttongroup  = a3buttongroup(buttons);

switch presentation
    case 'monaural'
        result=readfile_replace([dir_templates 'a3testacstimuli.xml'],experiment);
        % result=readfile_replace('a3testacstimuli-TF-20131029.xml',experiment);
    case 'diotic'
        result=readfile_replace([dir_templates 'a3testacstimuli-diotic.xml'],experiment);
end

fid=fopen([path savefilename '.apx302'],'w');
if (fid==-1)
    error(['Can''t open file ' savefilename ]);
end
fwrite(fid, result);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
