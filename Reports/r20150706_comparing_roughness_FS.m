function r20150706_comparing_roughness_FS
% function r20150706_comparing_roughness_FS
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%       See also r20141107_roughness.m, r20141126_roughness_validation.m
% 
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 05/07/2015
% Last update on: 05/07/2015 % Update this date manually
% Last use on   : 05/07/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

% % Roughness:
% outs1 = r20141107_roughness;
% outs3 = r20141126_roughness_validation;
% outs2 = r20150629_roughness_test;

options.bCreate = 0; % Set to 1 if you want to generate the test signals
options.bDoExp0 = 1;
options.bDoExp3 = 0;
options.bDoExp6 = 0;
% r20141126_roughness_validation(options);

% Fluctuation strength:
r20141126_fluctuation(options);

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
