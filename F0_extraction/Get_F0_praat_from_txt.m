function [t f0] = Get_F0_praat_from_txt(filename, info)
% function [t f0] = Get_F0_praat_from_txt(filename, info)
%
% Reads Fundamental frequency F0 from text file generated in Praat:
%   In Praat:
%       1. Open wav file
%       2. Praat objects window: press View and edit
%       3. On plots, select the whole file
%       4. Go to: Pitch -> Pitch listing
%       5. On the same window select Save as... 
%
% Desirable: to set the Pitch Range from 75 to 300 Hz (Pitch Settings...)
%
% % Example:
%   [t f0] = Get_F0_praat_from_txt('/home/alejandro/Documenten/MATLAB/MATLAB_svn/new_audio/wdz6.txt');
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, 2014
% Created in    : 2014
% Last update on: 19/05/2014
% Last use on   : 19/05/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    info = [];
end

info = Ensure_field(info, 'F0max', 1400);
info = Ensure_field(info, 'bPrint', 0);
info = Ensure_field(info, 'time2compensate', 0);

F0max = info.F0max;
bPrint = info.bPrint;
time2compensate = info.time2compensate;

if nargin < 3
    bPrint = 0;
end

literal         = '--undefined--'; % Character when no F0 is detected

count           = 1;

fid             = fopen(filename);
y               = 0;
tline           = fgetl(fid);

idx_t = 1;

while ischar(tline)&&length(tline)>2
    matches             = strfind(tline, literal);
    if isempty(matches);
        match_title     = strfind(tline, 'F0_Hz');
        if isempty(match_title) % then it is not a header
            tlinesplit = strsplit(tline,' ');
            try
                if length(tlinesplit{1}) == 0
                    idx_t = 2; % Then first character is null
                end
                t(count)    = str2num(tlinesplit{idx_t}); 
            catch
                t(count)    = str2num(tline(1:end-11)); % For ensuring compatibility with older versions of this script
            end
            try
                f0(count)   = str2num(tlinesplit{idx_t+1}); 
            catch
                f0(count)   = str2num(tline(end-10:end)); % For ensuring compatibility with older versions of this script
            end
            count       = count+1;
            idx_t = 1;
        end
    else
        t(count) = str2num(tline(1:end-14));
        f0(count) = NaN;
        count = count+1;
    end
         
    matches = strfind(tline, literal);
    num = length(matches);
    if num > 0
        y = y + num;
        if bPrint == 1
            fprintf(1,'%d:%s\n',num,tline);
        end
    end
    tline = fgetl(fid); % Get next line
end
fclose(fid);

idx = find(f0>F0max);
f0(idx) = NaN;
if length(idx) > 0
    disp([mfilename '.m: ' num2str(length(idx)) ' F0 values above F0max. The returned values for them are NaN'])
    pause(1)
end

t = t + time2compensate;
idx2delete      = find(t < 0);
t(idx2delete)   = [];
f0(idx2delete)  = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end