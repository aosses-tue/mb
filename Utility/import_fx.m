function [t F0] = import_fx(filename, info) 
% function [t F0] = import_fx(filename, info) 
%
%   info.startRow   - line number in txt-file where F0 starts being read, 
%                     default = 15, if *.fx files are loaded. 
%   info.time2compensate - 
%   filename = *.fx file, format according to Paul Bagshaw's database
%   t   - time in seconds
%   F0  - Fundamental frequency
%
% For an example, run script: F0mod_validation_PB_database.m
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2013-2014
% Created on: 2013
% Last update: 27/5/2014 % Update this date manually
% Last used: 27/5/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    info = [];
end
   
info = Ensure_field(info,'startRow', 15); % 15 in case of Paul Bagshaw's frequency contours
info = Ensure_field(info,'time2compensate',0);

startRow        = info.startRow;
time2compensate = info.time2compensate;

%% Initialize variables.
delimiter = {'\t',' '};
if nargin<=2
    startRow = 15;
end

if nargin == 0
    filename  = '~/Documenten/LaTeX_Docs/proefschrift-2014/Speech-material/fda_eval/sb/sb002.fx';
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec);

%% Close the text file.
fclose(fileID);

if length(dataArray) == 3
    dataArray{1,3}=[];
end

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end

raw_time_col = dataArray(:,1);
raw_F0_col   = dataArray(:,2);
 
R = cellfun(@(x) strcmp(x,'='), raw_time_col,'UniformOutput',false); % Find non-numeric cells
idx = find(cell2mat(R)==1);

t = str2double(raw_time_col{1,1});
t(idx) = NaN;

for k = 1:length(idx)
    try
        Delta = t( idx(k)-1 ) - t( idx(k)-2 );
        t( idx(k) ) = t( idx(k)-1 )+Delta;
    end
end
F0 = str2double(raw_F0_col{1,1}); % Replace non-numeric cells
F0(idx) = NaN;

t = t - time2compensate;
idx2delete      = find(t < 0);
t(idx2delete)   = [];
F0(idx2delete)  = [];

if max(t) > 500 % then I will assume time is in millisecons
    t = t/1000; % Conversion from ms to s
    txt_time = ' milliseconds';
else
    txt_time = ' seconds';
end

if time2compensate > 0
    disp([mfilename '.m: Compensating time in ' Num2Str(time2compensate) txt_time])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end