function misc = Start_TUe
% function misc = Start_TUe
%
% 1. Summary: 
%       - Initialises AMS Toolbox
%       - Gives several directories coded in the struct 'misc'
%
% 2. Additional info:
%   Tested cross-platform: Yes
%
% 3. Use:
%       If working connected to TUe server either physically or remotely, 
%       set info.bOnline to 1 otherwise set it to 0
% 
% 4. Stand-alone example:
%       Start_TUe;
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 01/05/2014
% Last update on: 28/07/2014
% Last use on   : 21/08/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc = Get_TUe_paths;

if nargout == 0
   
    Add_paths(misc);
    
    if isfield(misc,'tb_Loudness_v12')
        misc_sub = Get_TUe_subpaths('tb_Loudness_v12');
        disp(['Adding paths under: tb_Loudness_v12'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_NMT')
        misc_sub = Get_TUe_subpaths('tb_NMT');
        disp(['Adding paths under: tb_NMT'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_Plot4papers_JU')
        misc_sub = Get_TUe_subpaths('tb_Plot4papers_JU');
        disp(['Adding paths under: tb_Plot4papers_JU'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_Psysound32')
        % Initialises PsySound 3:
        psysound3(0);
    end
    
    amtstart;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename])