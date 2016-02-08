function txtout = Get_text_from_txt(file)
% function txtout = Get_text_from_txt(file)
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 08/02/2016
% Last update on: 08/02/2016 
% Last use on   : 08/02/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    file = 'D:\MATLAB_shared\Psychoacoustics\FluctuationStrength\Literature\Data\FM-fm.csv';
end

fid             = fopen(file);
tline           = fgetl(fid);

i = 1;
while ischar(tline)
    txtout{i} = tline;
    tline = fgetl(fid);
    i = i+1;
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
