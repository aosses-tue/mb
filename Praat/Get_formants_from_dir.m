function [t F0] = Get_formants_from_dir(folder, params)
% function [t F0] = Get_formants_from_dir(folder, params)
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       info.startRow = 2; 
%       [t formants] = Get_formants_from_dir('D:\MATLAB\Output\tmp-VoD_MIRtoolbox\modus-4_v3-2filt.txt',info);
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 12/08/2014
% Last update on: 12/08/2014 % Update this date manually
% Last use on   : 12/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    folder = 'D:\MATLAB\Output\tmp-VoD_MIRtoolbox\';
end

bDiary = 1;
Diary(mfilename,bDiary,folder);

if nargin < 2
    params = [];
end

params = Ensure_field(params,'timestep'    , 0.01); % positive timestep 0.01
params = Ensure_field(params,'nformants'   ,    5); % positive nformants 5
params = Ensure_field(params,'maxformant'  , 5500); % positive maxformant 5500
params = Ensure_field(params,'windowlength',0.025); % positive windowlength 0.025
params = Ensure_field(params,'dynamicrange',   20); % positive dynamic range 20

local_praat     = Get_TUe_paths('praat');
local_praat_sc  = Get_TUe_paths('praat_scripts');
script          = [local_praat_sc 'Get_formants_from_dir.praat'];

if isunix
    
else
	 command4system = [local_praat ' ' script ' ' folder ' ' ...
                                              '"' num2str(params.timestep)     '" ' ...
                                              '"' num2str(params.nformants)    '" ' ...
                                              '"' num2str(params.maxformant)   '" ' ...
                                              '"' num2str(params.windowlength) '" ' ...
                                              '"' num2str(params.dynamicrange) '"'];
end

disp([mfilename '.m: ' command4system])
[s r] = system( command4system );

if s == -1
    disp([mfilename '.m: problem running Praat, please check that praat is installed and the scripts do exist...'])
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
