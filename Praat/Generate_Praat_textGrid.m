function Generate_Praat_textGrid(filename,tlims)
% function Generate_Praat_textGrid(filename,tlims)
%
% 1. Description:
%       Generates textGrid to display in a nice way audio files in Praat
%       tlims specifies the limits of each audio segments.
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       paths = Get_TUe_paths('db_fastl2007');
%       filename = [paths 'track_38.wav'];
%       tlims = [1:3:41]; % I know in advance that the audio file is 41-seconds long
%       Generate_Praat_textGrid(filename,tlims);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 13/08/2014
% Last update on: 13/08/2014 
% Last use on   : 25/11/2014 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    filename = [Get_TUe_data_paths('db_fastl2007') 'track_38.wav'];
               % 'D:\Documenten-TUe\10-Referenties\02-Mijn-boeken\Fastl2007-psychoacoustics\Sound\
end

path_praat  = Get_TUe_paths('praat_scripts');
path_out    = Get_TUe_paths('outputs');
inputfile   = [path_praat 'template_textgrid.praat'];
try
    tmp = strsplit(filename,delim);
    inputwavfile = [path_out tmp{end}];
    outputfile = [Delete_extension( inputwavfile,'wav') '.TextGrid'];
catch
    [xx, tmp] = fileparts(filename);
    inputwavfile = [path_out tmp]; 
    outputfile = [inputwavfile '.TextGrid'];
end

[x fs] = Wavread(filename);
t = ( 1:length(x) )/fs;

if nargin < 2
    tlims = [1:3:max(t)];
end

if tlims(1) ~= 0
    tlims = [0 tlims];
end

if tlims(end) ~= max(t);
    tlims = [tlims max(t)];
end

toreplace = Get_date;
toreplace.wavfilename = Replace_character(inputwavfile,delim,'/');
toreplace.mfilename = mfilename;
toreplace.tmax = max(t); 
nframes = length(tlims)-1;
toreplace.line1 = sprintf('intervals: size = %i\n',nframes);

for i = 1:1:length(tlims)-1
    toreplace.line1 = sprintf('%sintervals [%i]:\nxmin = %.15f\nxmax = %.15f\ntext = "%s"\n',toreplace.line1,i,tlims(i),tlims(i+1),['s' num2str(i)]);
end

output = readfile_replace(inputfile,toreplace);

fid=fopen(outputfile, 'w');
fwrite(fid, output);
fclose(fid);

disp(['m-file: ' outputfile ' successfully created'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
