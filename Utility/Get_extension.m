function ext = Get_extension(text)
% function ext = Get_extension(text)
%
% Get_extension('filename.txt');
% ext will be 'txt'
%
% Programmed by Alejandro Osses, ExpORL, KULeuven, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ext = '';

temp = strsplit(text,'.');

if length(temp) >= 2
    ext = temp{end};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end