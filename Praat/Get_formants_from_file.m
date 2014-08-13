function [t F0] = Get_formants_from_file(filename, info)
% function [t F0] = Get_formants_from_file(filename, info)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       info.startRow = 2; 
%       [t formants] = Get_formants_from_file('D:\MATLAB\Output\tmp-VoD_MIRtoolbox\modus-4_v3-2filt.txt',info);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 11/08/2014
% Last update on: 11/08/2014 % Update this date manually
% Last use on   : 12/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % function [t F0] = import_fx(filename, info) 
        % % function [t F0] = import_fx(filename, info) 
        % %
        % %   info.startRow   - line number in txt-file where F0 starts being read, 
        % %                     default = 15, if *.fx files are loaded. 
        % %   info.time2compensate - 
        % %   filename = *.fx file, format according to Paul Bagshaw's database
        % %   t   - time in seconds
        % %   F0  - Fundamental frequency
        % %
        % % For an example, run script: F0mod_validation_PB_database.m
        % %
        % % Programmed by Alejandro Osses, ExpORL, KU Leuven, Belgium, 2013-2014
        % % Created on: 2013
        % % Last update: 27/5/2014 % Update this date manually
        % % Last used: 27/5/2014 % Update this date manually
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    info = [];
end
   
info = Ensure_field(info,'startRow', 2); % 15 in case of Paul Bagshaw's frequency contours
info = Ensure_field(info,'time2compensate',0);

startRow        = info.startRow;
time2compensate = info.time2compensate;

%% Initialize variables.
delimiter = {'\t',' '};
if nargin<=2
    startRow = 2;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
% formatSpec = '%s%s%[^\n\r]';
formatSpec = '%s%s%s%s%s%s%[^\n\r]';

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
try
    raw_F2_col = dataArray(:,3);
    raw_F3_col = dataArray(:,4);
    raw_F4_col = dataArray(:,5);
    raw_F5_col = dataArray(:,6);
end
 
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
try
    for i = 2:5
        Exp1 = ['F0 = [F0 str2double(raw_F' num2str(i) '_col{1,1})];'];
        eval(Exp1);
    end
end
F0(idx) = NaN;
F0 = Delete_NaN_columns(F0);

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

if nargout == 0
    plot(t,F0,'o');
    xlabel('Time')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
