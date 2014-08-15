function y = Generate_Fastl2007_textGrid
% function y = Generate_Fastl2007_textGrid
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 14/08/2014
% Last update on: 14/08/2014 % Update this date manually
% Last use on   : 14/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figures2generate = {    9, 3,'track_34'; ...
                       10, 1,'track_35'};
                   
for i = 1:size(figures2generate,1)
    tmpinfo = Check_Fastl2007(figures2generate{i,1}, figures2generate{i,2});
    Generate_Praat_textGrid([Get_TUe_paths('db_fastl2007') figures2generate{i,3} '.wav'],tmpinfo.timesep);
    
    % move to Fastl2007's directory
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
