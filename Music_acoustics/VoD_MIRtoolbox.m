function paths = VoD_MIRtoolbox(bDoEliminateClick)
% function paths = MIRtoolbox(bDoEliminateClick)
%
% 1. Description:
%       Applies MIT tb (different tools). For a demo on this, come back to 
%       any control version before (including) 06/08/2014.
% 
%       fs required if 'filename' is numeric.
%       MIR toolbox: interesting but discarded from our analysis according
%       to meeting with Armin on 23/07/2014 (no psychoacoustics). However
%       there are some very interesting tools.
% 
% 2. Additional info:
%       Tested cross-platform: No
%
% 3.1 Stand-alone example:
%       % Make sure you have not run PsySound3
%       % Make sure MIRToolbox v. 1.5 has been added to path
%       VoD_MIRtoolbox;
%
% % Audio excerpt into double:
%       xx = get(y,'Data');
%       xx = xx{1};
%       xx = xx{1};
%       
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 21/07/2014
% Last update on: 07/08/2014 % Update this date manually
% Last use on   : 07/08/2014 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    bDoEliminateClick = 1;
end

% if nargin == 0
    
if bDoEliminateClick == 1
    close all
end
m = Get_TUe_subpaths('db_voice_of_dragon');
rootfolder = m.dir_calibrated_m;
filenames = Get_filenames(rootfolder,'*2filt.wav'); % near field
                                                    % Ac.mode = 2, measured signal

disp('In case is not running, please remove MIR tb from PsySound and add MIR toolbox v15')
% disp('press any button to continue...')
pause(2)
% addpath('D:\MATLAB_git\tb_MIR_v15\MIRToolbox\')

% end

misc    = Get_VoD_params(1);
nPeriods = 10;

outputdir = [Get_TUe_paths('outputs') 'tmp-' mfilename delim];
if bDoEliminateClick
    Mkdir(outputdir);
end

for mode_idx = 1:length(filenames)
    
    ac_mode     = mode_idx + 1;
    filename    = [rootfolder filenames{mode_idx}];
    outfilename{mode_idx} = [outputdir filenames{mode_idx}];
    ti          = misc.ti_measured(mode_idx);
    tf          = ti + nPeriods*misc.Tmodel(mode_idx);
    t_before_after = 0.025*misc.Tmodel(mode_idx); % time before and after theoretical click time, total of 5%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if bDoEliminateClick % We eliminate click and generate temporal wav-files
        
        figure;
        miraudio(filename, 'Extract',ti,tf)
        y = miraudio(filename, 'Extract',ti,tf); % extracting from second 1 to 2
        ytmp    = Get_MIR_param(y,'Data');
        fs      = Get_MIR_param(y,'Sampling');
        
        [z tinfo] = Zerocross(y);

        for i = nPeriods:-1:1

            t_per = misc.t_coil(ac_mode-1,i);
            idx1 = max( find(tinfo.tup < t_per-t_before_after) );
            tup = tinfo.tup(idx1);
            idx1 = find(tup == tinfo.t); % relative to time vector
            
            idx2 = min( find(tinfo.tdown > t_per+t_before_after) );
            tdown = tinfo.tdown(idx2);
            idx2 = find(tdown == tinfo.t);
            
            yclean = ytmp([1:idx1,idx2:end]);
            
            fprintf('To delete from %.4f to %.4f [s], click at = %.4f [s]\n',tup,tdown,t_per);
            
        end
        
        Wavwrite(yclean,fs,outfilename{mode_idx});
        
    end
        
end
    
paths = outfilename;

% Plot_vline_periodic(misc.     )

% % Roughness:
% mirroughness(y)
% % mircentroid(a,'Frame')
%  
% % Fluctuation
% mirfluctuation(y)
% mirfluctuation(y,'Summary')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end