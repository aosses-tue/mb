function Save(filename,variable)
% function Save(filename,variable)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 21/10/2014
% Last update on: 21/10/2014 % Update this date manually
% Last use on   : 21/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = fileparts(filename);

if strcmp(str,'')
    try
        str = Get_TUe_paths('outputs');
    catch
        str = [pwd delim];
    end
    filename = [str filename];
end

save(filename,variable);
disp(['Variable saved as: ' filename '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
