function status = Mkdir(folder)
% function status = Mkdir(folder)
%
% 1. Description:
%       Creates a folder only in case it does not exist previously
%       'status' is 1 if the directory was sucessfully created
%       'status' is 0 if the directory already exists.
%
% Programmed by Alejandro Osses, ExpORL, KU Leuven, 2014
% Last update on: 12/08/2014
% Last use on   : 16/03/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bCreateFolder = ~isdir(folder);

if bCreateFolder == 1
    mkdir(folder);
    disp(['Directory ' folder 'created successfully...'])
    status = 1;
else
    if nargout == 1
        warning('Directory not created, it already exists, you might be overwriting files');
    end
    status = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
