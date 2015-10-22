function misc = Start_TUe
% function misc = Start_TUe
%
% 1. Description: 
%       - Initialises AMT Toolbox
%       - Gives several directories coded in the struct 'misc'
%
% 2. Use:
%       If working connected to TUe server either physically or remotely, 
%       set info.bOnline to 1 otherwise set it to 0
% 
% 3. Stand-alone example:
%       Start_TUe;
%
% 4. Additional info:
%   Tested cross-platform: Yes
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 01/05/2014
% Last update on: 14/08/2015
% Last use on   : 21/10/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

misc = Get_TUe_paths;

if nargout == 0
   
    Add_paths(misc);
    
    if isfield(misc,'Psychoacoustics')
        misc_sub = Get_TUe_subpaths('Psychoacoustics');
        disp(['Adding paths under: Psychoacoustics'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_AM_AddOns')
        misc_sub = Get_TUe_subpaths('tb_AM_AddOns');
        disp(['Adding paths under: tb_AM_AddOns'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_Loudness_v12')
        misc_sub = Get_TUe_subpaths('tb_Loudness_v12');
        disp(['Adding paths under: tb_Loudness_v12'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_Misc')
        misc_sub = Get_TUe_subpaths('tb_Misc');
        disp(['Adding paths under: tb_Misc'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_NMT')
        misc_sub = Get_TUe_subpaths('tb_NMT');
        disp(['Adding paths under: tb_NMT'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_NMTAddOns')
        misc_sub = Get_TUe_subpaths('tb_NMTAddOns');
        disp(['Adding paths under: tb_NMTAddOns'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'tb_Plot4papers_JU')
        misc_sub = Get_TUe_subpaths('tb_Plot4papers_JU');
        disp(['Adding paths under: tb_Plot4papers_JU'])
        Add_paths(misc_sub);
    end
    
    if isfield(misc,'Reports_KUL')
        misc_sub = Get_TUe_subpaths('Reports_KUL');
        disp(['Adding paths under: Reports_KUL'])
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