function [m,path,dBFS,noisefile,info]=es_speech_material(name, reread)
% function [m,path,dBFS,noisefile,info]=es_speech_material(name, reread)
%
% m is a cell with 3 dimensions: set, list, word/sentence
% dBFS is the average signal level on disk
% Noisefile is the wav file containing the noise to be used with the
% material (usually LTASS)
%
%   info.speechmaterial_dirlists -  folder where all the speech material 
%                                   (wav files are)
%
% 
% Example:
%       name = 'Matrix';
%       reread = 0;
%       [m,path,dBFS,noisefile,info]=es_speech_material(name, reread);  
%
% Created by    : Alejandro Osses, TUe 2015
% Created on    : 14/06/2015
% Original file : nl_speech_material.m
% Last use on: 14/06/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = [];

if nargin==0
    return;
end

if (nargin<2)
    reread=0;
end

mpath=fileparts(mfilename('fullpath'));
cachefile=[mpath '/es_speech_material.mat'];
if (~reread && exist(cachefile, 'file'))
    load(cachefile);
    
    if (isfield(cache, name))
        m=cache.(name).m;
        path=cache.(name).path;
        dBFS=cache.(name).dBFS;
        noisefile=cache.(name).noisefile;
        info=cache.(name).info;
        return;
    end
end

noisefile = '';
switch name
    case 'Matrix'
        [m,path,info]=esmatrix_l;
        dBFS = -27; % you can confirm with rmsdb('noisevalue')
        noisefile = 'EsMatrixnoise_ltass.wav';
    otherwise 
        error('Invalid material');
end

% Save cache
if (exist(cachefile, 'file'))
    load(cachefile);
else
    cache=struct;
end
cache.(name).m=m;
cache.(name).path=path;
cache.(name).dBFS=dBFS;
cache.(name).noisefile=noisefile;
cache.(name).info=info;
save(cachefile, 'cache');

% cell(set,list,sentence)

% s.file
% s.text
% s.keywords

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,esmatrixpath,info]=esmatrix_l
% vlaamse matrix, basislijsten
%
% output example: 
%   - m{1,10,12}
%       text:       'The milk was by the front door'
%       keywords:   '2 6 7'
%       file:       'BKB10 S12.wav'
%
% esmatrixpath:
%   - '/home/alejandro/x/speechmaterials/dutch/Matrix/'
%
% Update also 'experiment_calibration_APEX_F0mod.m' to load the proper noise
% file if necessary

nLists      = 26;
nSentences  = 10;
txt_file    = 'EsMatrixBasisList.txt';
m = cell(1,nLists,nSentences);
try
    info.path = Get_TUe_paths('db_speechmaterials');
catch
    pa = getpaths; % Tom's configuration
    info.path = pa.speechmaterials;
end

info.uri = 'spanish/Matrix/';
info.speechmaterial_dirlists = '';

% vlmatrixpath=[pa.local_speechmaterials info.uri]; % x-Drive, ExpORL
esmatrixpath=[info.path info.uri]; % x-Drive, ExpORL

fid=fopen([esmatrixpath txt_file]);

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if (~isempty(regexp(tline, '^\d+$', 'once')))
        listnr=str2num(tline);
        sentnr=0;
        continue;
    end
    
    %                      (not digit) followed by (digits)
    tokens=regexp(tline, '^(.*\D) ([0-9][0-9 ]+)', 'tokens');
    if (isempty(tokens))
        tokens{1}{1} = tline;
        tokens{1}{2} = '1 1';
        warning(['Could not parse line ' tline ', probably no keywords were found']);
    end
    sentnr=sentnr+1;
    s=struct;
    s.text = tokens{1}{1};
    s.keywords = tokens{1}{2};
    filename = Get_EsMatrix_wav_filename(s.text);
    s.file = [filename '.wav'];
    m{1,listnr,sentnr}=s;
end

fclose(fid);

disp('');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function news=replstr_l(s, from_to)

news=[];
for n=1:length(s)
    found=0;
    for i=1:size(from_to,1)
        if (strncmp(s(n:end), from_to{i,1}, length(from_to{i,1})))
            news=[news from_to{i,2}];
            n=n+length(from_to{i,1})-1;
            found=1;
            break;
        end
    end
    if (~found)
        error('Phoneme not found');
    end
end
