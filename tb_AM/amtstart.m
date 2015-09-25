function amtstart(varargin)
% function amtstart(varargin)
% 
% 1. Description:
%       Start the Auditory Modeling Toolbox
% 
%   Usage:  amtstart;
%           amtstart(flags);
% 
%   AMTSTART starts the Auditory Modeling Toolbox. This command must be
%   run before using any of the function in the toolbox.
%
%   The AMT depends on the Linear Time Frequency Analysis Toolbox (LTFAT). 
%   You must first download LTFAT from
%   http://ltfat.sourceforge.net/ and unpack the downloaded file. 
%   In the AMT, there is a pre-prepared directory /thirdparty/ltfat
%   where the LTFAT can be stored. Alternatively, set the path to your
%   LTFAT installation to the search path of Matlab/Octave.
%
%   In order to run all the models from AMT, you will need:
%   
%   1) run amtmex and compile successfully
%   2) SOFA API from http://sourceforge.net/projects/sofacoustics for Matlab (in e.g. amtbase/thirdparty/SOFA)
%   3) Python >2.6 is required with numpy and scipi packages. On Linux, use sudo apt-get install python-scipy python-numpy
%   4) run make (Linux) or make.bat (Windows) in amtbase/src/verhulst
%   5) Optimization Toolbox for Matlab
%   6) Data in amtbase/signals/ and amtbase/hrtf/ depending on the model
%
%   See also:  amthelp
%
%   Url: http://amtoolbox.sourceforge.net/doc/amtstart.php

% Copyright (C) 2009-2014 Peter L. Søndergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%   AUTHOR : Peter L. Soendergaard, Piotr Majdak 
% Last update: 25/09/2015 by AO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bp=amtbasepath; % bp - basepath

% Search for LTAFT package

if ~exist('ltfatstart','file')
  ltfatpath=fullfile(bp,'thirdparty','ltfat');
  if exist(ltfatpath,'dir')
    addpath(ltfatpath);
  end
end

% Load the version number
[FID, MSG] = fopen ([bp,'amtoolbox_version'],'r');
if FID == -1
    error(MSG);
else
    amt_version = fgetl (FID);
    fclose(FID);
end

% Check if 'silent' present in the flags
silent=0;
if exist('OCTAVE_VERSION','builtin'), args=argv; else args=varargin; end
for ii=1:numel(args)
    s=lower(args{ii});
    if strcmp(s,'silent') || strcmp(s,'-q')
        silent=1;
    end;
end;

if ~silent
    disp('  ');
    disp(['AMT version ',amt_version,'. (C) Piotr Majdak and Peter L. Soendergaard.']);
    disp('  ');
    disp('Starting toolboxes...');
end;
 
%% LTFAT package

if exist('ltfatstart','file')
  if silent, ltfatstart(0); else ltfatstart; end;
else
  error(['LTFAT package could not be found. Unable to continue.' 10 ...
        'Download LTFAT from http://ltfat.sourceforge.net ' 10 ...
        'and copy to amtoolbox/thirdparty/ltfat.']); 
end

% Check for the correct version. 
s=ltfathelp('version'); 
s_r='1.0.9'; % set the required version
warning('Update LTFAT version')
% % For AMT version 0.9.7
% s_r='2.0.0'; % set the required version
v=sscanf(s,'%d.%d.%d'); v(4)=0;
v_r=sscanf(s_r,'%d.%d.%d');
if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) ),
    error(['You need LTFAT >= ' s_r ' to work with AMT. ' ...
      'Please update your package from http://ltfat.sourceforge.net ']);
end  
     
%% SOFA package

% Search for SOFA package
if ~exist('SOFAstart','file')
    sofapath=fullfile(bp,'thirdparty','SOFA','API_MO');
    if exist(sofapath,'dir')
        addpath(sofapath);
    end
end

% Start SOFA
if exist('SOFAstart','file')
    SOFAdbPath(fullfile(bp,'hrtf'));
    SOFAdbURL('http://www.sofacoustics.org/data'); % This is a default path and will be overwritten later
    if silent, SOFAstart('silent'); else SOFAstart('short'); end
    warning('off','SOFA:upgrade');	% disable warning on upgrading older SOFA files
	warning('off','SOFA:load'); % disable warnings on loading SOFA files
else
  disp(['SOFA package could not be found. Continue without SOFA support.']);
  disp(['For SOFA support please download the package ' ...
        'from http://sofacoustics.sourceforge.net ' ...
        'and copy to amtoolbox/thirdparty/SOFA.']); 
    if ~silent,
        disp('SOFA package could not be found. Continue without SOFA support.');
        disp(['For SOFA support please download the package ' ...
            'from http://sofacoustics.sourceforge.net ' ...
            'and copy to amtoolbox/thirdparty/SOFA.']); 
    end
end

%% SFS package

% Search for the package
bp=which('amtstart');
bp=bp(1:end-11);
if ~exist('SFS_start','file')
  sfspath=fullfile(bp,'thirdparty','SFS');
  if exist(sfspath,'dir')
    addpath(sfspath);
  end
end

% Start 
disp('*** Starting SFS ***');
if exist('SFS_start','file')
  SFS_start;
  s=SFS_version; s_r='0.2.4'; % set the required version
  disp(['SFS, version ' s]);
  v=sscanf(s,'%d.%d.%d'); v(4)=0;
  v_r=sscanf(s_r,'%d.%d.%d');
  if ~(v(1)>v_r(1) || (v(1)>=v_r(1) && v(2)>v_r(2)) || (v(1)>=v_r(1) && v(2)>=v_r(2) && v(3)>=v_r(3)) ),
      error(['You need SFS >= ' s_r ' to work with AMT. ' ...
        'Please update your package from https://github.com/sfstoolbox/sfs ']);
  end  
else
  disp(['SFS package could not be found. Continue without SFS support.']);
  disp(['For SFS support please download the package ' ...
        'from https://github.com/sfstoolbox/sfs ' ...
        'and copy to amtoolbox/thirdparty/SFS.']); 
end

%% Start AMT
disp('*** Starting AMT ***');  
% --- general settings ---
% Print the banner at startup?
printbanner=1;

% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
bp=which('amtstart');
bp=bp(1:end-11);
if exist('addpath','var')>0
  addpath(bp);
else
  path(path,bp);
end
bp=[bp,filesep];

% Load the version number
[FID, MSG] = fopen ([bp,'amtoolbox_version'],'r');
if FID == -1
    error(MSG);
else
    amt_version = fgetl (FID);
    fclose(FID);
end

% -----------  install the modules -----------------

modules={};
nplug=0;

% List all files in base directory
d=dir(bp);

for ii=1:length(d)
  if d(ii).isdir
    if ~(d(ii).name(1)=='.')

      name=d(ii).name;
      
      % The file is a directory and it does not start with '.' This could
      % be a module      
      if exist([bp,name,filesep,name,'init.m'],'file')
	% Set 'status' to zero if the module forgets to define it.
	status=0;
	module_version=amt_version;
        addpath([bp,name]);

	eval([name,'init']);
        if status>0
          if status==1
            nplug=nplug+1;
            modules{nplug}.name=name;
            modules{nplug}.version=module_version;
          end;
	else
	  rmpath([bp,name]);
	end;
      end;	

    end;
  end;
end;

% Check if Octave was called using 'silent'
%if isoctave
%  args=argv;
%  for ii=1:numel(args)
%    s=lower(args{ii});
%    if strcmp(s,'--silent') || strcmp(s,'-q')
%      printbanner=0;
%    end;
%  end;
%end;

if printbanner
  disp(['AMT version ',amt_version,'. (C) Peter L. Søndergaard and Piotr Majdak. For help, please type "amthelp".'])
end;



%% ---------- load information into ltfathelp ------------

% As comp is now in the path, we can call ltfatarghelper
ltfatsetdefaults('amthelp','versiondata',amt_version,...
                 'modulesdata',modules);

