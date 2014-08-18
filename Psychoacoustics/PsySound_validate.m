function PsySound_validate
% function PsySound_validate
%
% 1. Description:
%
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 17/08/2014
% Last update on: 17/08/2014 % Update this date manually
% Last use on   : 17/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 1;
Diary(mfilename,bDiary);

% use: Generate_reference_sounds to generate reference sounds

dir_where_ref = Get_TUe_paths('outputs');

ref_loud  = [dir_where_ref 'ref_loud.wav'];
ref_sharp = [dir_where_ref 'ref_sharp.wav'];
ref_fluct = [dir_where_ref 'ref_fluct.wav'];
ref_rough = [dir_where_ref 'ref_rough.wav'];

% 1. Loudness
options.nAnalyser = 12; % DLM Chalupper
options.calfile = '~/Documenten/MATLAB/outputs/track_03.wav';
options.callevel = 60;
[tmp_h tmp_ha]  = PsySoundCL(ref_loud,options);

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
