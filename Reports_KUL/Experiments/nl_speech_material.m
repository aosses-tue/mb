function [m,path,dBFS,noisefile,info]=nl_speech_material(name, reread)
% function [m,path,dBFS,noisefile,info]=nl_speech_material(name, reread)
%
% m is a cell with 3 dimensions: set, list, word/sentence
% dBFS is the average signal level on disk
% Noisefile is the wav file containing the noise to be used with the
% material (usually LTASS)
%
%   info.speechmaterial_dirlists -  folder where all the speech material 
%                                   (wav files are)
%
% Adapted from oz_speech_materials by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info = [];

if nargin==0
    return;
end

if (nargin<2)
    reread=0;
end

mpath=fileparts(mfilename('fullpath'));
cachefile=[mpath '/nl_speech_material.mat'];
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
    case 'BKB'
        [m,path,info]=bkb_l;
        dBFS = -27;     % avgsignalrms: -27.3652
    case 'CNC_male'
        [m,path]=cnc_l;
        path = [path 'ALLRD/'];
        dBFS = -27;      % avgsignalrms: -27.1461  
    case 'CNC_female'
        [m,path]=cnc_l;
        path = [path 'ALLKG/'];
        dBFS = -27;      % avgsignalrms: -27.2525
    case 'Cuny'
        [m,path]=cuny_l;
        dBFS = -25;       % Based on noise track
    case 'audiobook'
        [m,path]=audiobook_l;
        dBFS = -27;       % Guesstimate
    case 'VCV'
        [m,path]=vcv_l;
        dBFS = -25;
    case 'LISTf'
        [m,path,info]=listf_l;
        dBFS = -25.2;
    case 'LISTm'
        [m,path,info]=listm_l;
        dBFS = -27;
    case 'Lilliput'
        dBFS = -28.4;
    case 'VlMatrix'
        [m,path,info]=vlmatrix_l;
        dBFS = -27; % you can confirm with rmsdb('noisevalue')
        noisefile = 'VlMatrixnoise_ltass.wav';
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
function [m,bkbpath,info]=bkb_l

% output example: 
%   - m{1,10,12}
%       text:       'The milk was by the front door'
%       keywords:   '2 6 7'
%       file:       'BKB10 S12.wav'
%
% bkbpath:
%   - '/home/alejandro/x/speechmaterials/australian/BKB/'

pa=getpaths;
info.path = pa.speechmaterials;
info.uri = 'australian/BKB/';
info.speechmaterial_dirlists = '';

%bkbpath=[pa.speechmaterials 'english/BKB/'];
bkbpath=[info.path info.uri]; % x-Drive, ExpORL

nLists      = 21;
nSentences  = 16;

m = cell(1,nLists,nSentences);

fid=fopen([bkbpath 'BKBList.txt']);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if (~isempty(regexp(tline, '^\d+$', 'once')))
        listnr=str2num(tline);
        sentnr=0;
        continue;
    end
    
    tokens=regexp(tline, '^(.*\D) ([0-9][0-9 ]+)', 'tokens');
    if (isempty(tokens))
        warning(['Could not parse line ' tline]);
    end
    sentnr=sentnr+1;
    s=struct;
    s.text = tokens{1}{1};
    s.keywords = tokens{1}{2};      % fixme
    s.file = sprintf('BKB%02d S%02d.wav', listnr, sentnr);
    m{1,listnr,sentnr}=s;
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,vlmatrixpath,info]=vlmatrix_l
% vlaamse matrix, basislijsten
%
% output example: 
%   - m{1,10,12}
%       text:       'The milk was by the front door'
%       keywords:   '2 6 7'
%       file:       'BKB10 S12.wav'
%
% vlmatrixpath:
%   - '/home/alejandro/x/speechmaterials/dutch/Matrix/'
%
% Update also 'experiment_calibration_APEX_F0mod.m' to load the proper noise
% file if necessary

nLists      = 26;
nSentences  = 10;
txt_file    = 'VlMatrixBasisList.txt';
m = cell(1,nLists,nSentences);
pa=getpaths;

info.path = pa.speechmaterials;
info.uri = 'dutch/Matrix/';
info.speechmaterial_dirlists = '';

vlmatrixpath=[pa.local_speechmaterials info.uri]; % x-Drive, ExpORL

fid=fopen([vlmatrixpath txt_file]);

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
    s.keywords = tokens{1}{2};      % fixme
    filename = Get_VlMatrix_wav_filename(s.text);
    s.file = [filename '.wav'];
    m{1,listnr,sentnr}=s;
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,listpath,info]=listf_l
% LIST-f
% Update also 'experiment_calibration_APEX_F0mod.m' to load the proper noise
% file if necessary

nLists      = 35;
nSentences  = 10;
txt_file    = 'LIST-fList.txt';
m           = cell(1,nLists,nSentences);
pa          = getpaths;

info.path   = pa.speechmaterials;
info.uri    = 'dutch/list/';
info.speechmaterial_dirlists = 'lijst';

listpath    = [pa.local_speechmaterials info.uri]; % x-Drive, ExpORL

fid         = fopen([listpath txt_file]);

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if (~isempty(regexp(tline, '^\d+$', 'once')))
        listnr=str2num(tline);
        sentnr=0;
        continue;
    end
    tokens=regexp(tline, '^(.*\D) ([0-9][0-9 ]+)', 'tokens');
    if (isempty(tokens))
        tokens{1}{1} = tline;
        tokens{1}{2} = '1 1';
        warning(['Could not parse line ' tline ', probably no keywords were found']);
    end
    sentnr=sentnr+1;
    s=struct;
    s.text = tokens{1}{1};
    s.keywords = tokens{1}{2};      % fixme
    
    filename = Get_LISTf_wav_filename(s.text,info);
    
    tmp_text = Replace_character(s.text,'ë','e');
    if ~strcmp(tmp_text,s.text)
        warning(['Text: ' s.text ', related to wav file ' filename ' contains a non compatible character with APEX3 experiments, it will be change'])
        s.text = tmp_text;
        pause(2)
    end
    
    s.file = [filename '.wav'];
    m{1,listnr,sentnr}=s;
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,listpath,info]=listm_l
% LIST-m
% Update also 'experiment_calibration_APEX_F0mod.m' to load the proper noise
% file if necessary

nLists      = 38;
nSentences  = 10;
txt_file    = 'LIST-mList.txt';
m           = cell(1,nLists,nSentences);
pa          = getpaths;

info.path   = pa.speechmaterials;
info.uri    = 'dutch/LISTman/'; % kijk maar: infilename=[speechpath p.speechmaterial_dirlists '/' sm.file];
info.speechmaterial_dirlists = ''; % Wav-files directly inside folder (without subfolders)

listpath    = [pa.local_speechmaterials info.uri]; % x-Drive, ExpORL

fid         = fopen([listpath txt_file]);

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if (~isempty(regexp(tline, '^\d+$', 'once')))
        listnr=str2num(tline);
        sentnr=0;
        continue;
    end
    tokens=regexp(tline, '^(.*\D) ([0-9][0-9 ]+)', 'tokens');
    if (isempty(tokens))
        tokens{1}{1} = tline;
        tokens{1}{2} = '1 1';
        warning(['Could not parse line ' tline ', probably no keywords were found']);
    end
    sentnr=sentnr+1;
    s=struct;
    s.text = tokens{1}{1};
    s.keywords = tokens{1}{2};      % fixme
    filename = Get_LISTm_wav_filename(s.text,info);
    
    tmp_text = Replace_character(s.text,'ë','e');
    if ~strcmp(tmp_text,s.text)
        warning(['Text: ' s.text ', related to wav file ' filename ' contains a non compatible character with APEX3 experiments, it will be change'])
        s.text = tmp_text;
        pause(2)
    end
   
    s.file = [filename '.wav'];
    disp(s.file);
    m{1,listnr,sentnr}=s;
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m,cncpath]=cnc_l
pa=getpaths;
cncpath=[pa.speechmaterials 'english/CNC/'];

m=cell(1,30,50);        % set, list, word

fid=fopen([cncpath 'Word2Phoneme.txt']);
w2p=struct;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end

    c=strsplit(tline, tb);
    c{1}=strrep(c{1}, '(', '_');
    c{1}=strrep(c{1}, ')', '_');
    c{1}=strrep(c{1}, '-', '_');
    c{1}=strrep(c{1}, ' ', '');
    c{3}=strrep(c{3}, ' ', '');
    w2p.(c{1}) = c{3};
end
fclose(fid);

fid=fopen([cncpath 'macquarie2CNC.txt']);
m2c=cell(10,2);       % set, list, word
count=1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if (isempty(tline))
        continue;
    end

    c=strsplit(tline, tb);
    m2c{count,1} = c{1};
    m2c{count,2} = c{2};
    count=count+1;
end
fclose(fid);


listnr=0;
fid=fopen([cncpath 'CNCLists.txt']);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if (isempty(tline))
        continue;
    end
    if (~isempty(regexp(tline, '^\d+$', 'once')))       
        listnr=str2num(tline);
        sentnr=0;
        continue;
    end
    
    sentnr=sentnr+1;
    
    list=rem(listnr,100);
    list=list + ((listnr-list)/100-1)*10;
    s=struct;
    s.text = tline;
    s.text=strrep(s.text, ' ', '');
    s.phoneme = upper(w2p.(s.text));
    s.cncphoneme = replstr_l(s.phoneme, m2c);
%     s.text = tokens{1}{1};
%     s.keywords = tokens{1}{2};      % fixme
    s.file = [s.phoneme '.WAV'];
    m{1,list, sentnr}=s;
    
    % Check if file exists
    if (~exist([cncpath '/ALLKG/' s.file], 'file'))
        error(['File ' s.file ' not found']);
    end
end

for s=1:size(m,1)
    for l=1:size(m,2)
        for q=1:size(m,3)
            if (isempty(m{s,l,q}))
                fprintf('CNC set %d, list %d, word %d missing\n', s, l, q);
            end
        end
    end
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,cunypath]=cuny_l


pa=getpaths;
cunypath=[pa.speechmaterials 'english/Cuny/'];

m=cell(1,170,12);     % list, sentence

fid=fopen([cunypath 'CunyList.txt']);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if (~isempty(regexp(tline, '^\d+$', 'once')))
        listnr=str2num(tline);
        sentnr=0;
        continue;
    end
    
    if (isempty(tline))
        continue;
    end
    sentnr=sentnr+1;
    s=struct;
    s.text = tline;

    subnr=floor((listnr-1)/10)*10;

    s.file = sprintf('CUNY%03d/L%03dS%02d.WAV', subnr ,listnr, sentnr);
    m{1,listnr,sentnr}=s;
    
    if (~exist([cunypath '/' s.file], 'file'))
        error(['File ' s.file ' not found']);
    end
end
fclose(fid);

for l=1:size(m,2)
    for s=1:size(m,3)
        if (isempty(m{1,l,s}))
            fprintf('Cuny list %d sentence %d missing\n', l, s);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,bookpath]=audiobook_l

pa=getpaths;
bookpath=[pa.speechmaterials 'english/audiobook/'];

m=cell(1,1,6);     % set, list, sentence

for i=1:size(m,1)
    for j=1:size(m,2)
        for k=1:size(m,3)
            s.file = sprintf('set%d/list%d/sentence%d.wav', i,j, k);
            m{i,j,k} = s;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m,bookpath]=vcv_l

pa=getpaths;
bookpath=[pa.speechmaterials 'english/VCV/'];

m=cell(1,2,16);     % set, list, sentence

for i=1:size(m,1)
    for j=1:size(m,2)
        switch j
            case 1
                prefix = 'male_a';
            case 2
                prefix = 'male_i';
        end
        
        files = dir([bookpath prefix '/*.wav']);
        
        for ifile=1:length(files)
            filename = files(ifile).name;
        
            s.file = sprintf('%s/%s', prefix, filename);
            s.text = filename(2:end-4);
            m{i,j,ifile} = s;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%