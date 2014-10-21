function r20141010_update_dau1996
% function r20141010_update_dau1996
%
% 1. Description:
%       Script used to generate plots included in report on 10/10/2014.
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3. Stand-alone example:
%       demo_dau1996b;
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 10/10/2014
% Last update on: 10/10/2014 % Update this date manually
% Last use on   : 20/10/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

demo_dau1996b; % as backed-up on 10/10/2014

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
