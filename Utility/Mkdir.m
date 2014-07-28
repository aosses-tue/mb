function Mkdir(folder)
% function Mkdir(folder)
%
% Creates a folder only in case it does not exist previously
%
% Programmed by Alejandro Osses, ExpORL, KULeuven, 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bCreateFolder = ~isdir(folder);

if bCreateFolder == 1
    mkdir(folder);
    disp(['Directory ' folder 'created successfully...'])
else
    warning('Directory not created, it already exists, you might be overwriting files');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end