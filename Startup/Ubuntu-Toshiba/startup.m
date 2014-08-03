% startup.m
% Location: '/home/alejandro/startup.m'
%
% Programmed by Alejandro Osses, TU/e 2014
% PC:           Toshiba (mi compu), Ubuntu
% Created on    : 18/05/2014
% Last update on: 04/08/2014
% Last use on   : 04/08/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

info.bOnline = 0; % default = 0
bTUe = 1;
bKUL = ~bTUe;

if bTUe
    if info.bOnline == 1
        % Last used: xx/xx/xxxx
        addpath('~/Documenten/Dropbox/TUe/MATLAB/Utility');
        addpath('~/Documenten/Dropbox/TUe/MATLAB/tb_NMT_4.31/Matlab/Processing/');
        disp('Using on-line set-up');
    else
        % Last used: 04/08/2014
        % addpath('~/Documenten/MATLAB/MATLAB_TUe/Utility');
        % addpath('~/Documenten/MATLAB/MATLAB_TUe/tb_NMT_4.31/Matlab/Processing/');
        addpath('~/Documenten/MATLAB/MATLAB_git/Utility');
        addpath('~/Documenten/MATLAB/MATLAB_git/tb_NMT_4.31/Matlab/Processing/');
        disp('Using local TU/e MATLAB files, change bTUe to 0 if you want to use the KUL set-up');
    end
    Start_TUe;
elseif bKUL
    % Last used: 28/06/2014
    addpath('~/Documenten/MATLAB/MATLAB_svn/Utility/');
    localpaths = Get_paths([],1); 
    disp('Using KU Leuven set-up, change bTUe to 1 to revert this option')
end

slCharacterEncoding('UTF-8');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([mfilename '.m: startup file for ''my Toshiba'' PC, using Ubuntu'])
