function Save(filename,variable,varargin)
% function Save(filename,variable,varargin)
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

if nargin == 2
    save(filename,variable);
else
    exp1 = ['save(filename,variable'];
    for i = 1:length(varargin)
        exp1 = [exp1 ',varargin{' num2str(i) '}'];
    end
    exp1 = [exp1 ');'];
        
    eval(exp1);
end
disp(['Variable saved as: ' filename '.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
