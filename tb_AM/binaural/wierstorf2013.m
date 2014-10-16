function [localization_error,perceived_direction,desired_direction,x,y,x0] = ...
        wierstorf2013(X,Y,phi,xs,src,L,varargin);
%WIERSTORF2013 estimate the localization within a WFS or stereo setup
%   Usage: [...] = wierstorf2013(X,Y,phi,xs,src,L,resolution,method,...);
%
%   Input parameters:
%       X           : range of the x-axis [xmin,xmax] in m (for a single x point
%                     you have to specify [x x]
%       Y           : range of the y-axis [ymin ymax] in m (for a single y point
%                     you have to specify [y y]
%       phi         : orientation of the listener in rad (0 is in the direction
%                     of the x-axis
%       xs          : position of the point source in m / direction of the 
%                     plane wave
%       src         : source type
%                       'ps' for a point source
%                       'pw' for a plane wave
%       L           : length/diameter of the loudspeaker array
%       resolution  : resolution x resolution is the number of points the
%                     localization should be estimated in the listening area.
%                     The points are evenly distributed along the axes. 
%       method      : reproduction setup
%                       'wfs' for wave field synthesis
%                       'setreo' for stereophony
%
%   Output parameters:
%       localization_error  : deviation from the desired direction, defined as
%                             perceived_direction - desired_direction / rad
%       perceived_direction : the direction of arrival the binaural model has
%                             estimated for our given setup / rad
%       desired_direction   : the desired direction of arrival indicated by the
%                             source position xs / rad
%       x                   : corresponding x-axis
%       y                   : corresponding y-axis
%       x0                  : position and directions of the loudspeakers in the
%                             form n x 6, where n is the number of loudspeakers
%
%   WIERSTORF2013(X,Y,phi,xs,src,'wfs',L,nls,array) calculates the localization
%   error for the defined wave field synthesis or stereophony setup. The
%   localization error is defined here as the difference between the perceived
%   direction as predicted by the dietz2011 binaural model and the desired
%   direction given by xs. The loudspeaker setup for the desired reproduction
%   method is simulated via HRTFs which are than convolved with white noise
%   which is fed into the binaural model.
%
%   The following parameters may be passed at the end of the line of
%   input arguments:
%
%     'resolution',resolution
%                      Resolution of the points in the listening
%                      area. Number of points is resoluation x resolution. If
%                      only one point in the listening area is given via single
%                      values for X and Y, the resolution is automatically set
%                      to 1.
%
%     'nls',nls        Number of loudspeaker of your WFS setup.
%                      Default value is 2.
%
%     'array',array    Array type to use, could be 'linear' or 'circle'.
%                      Default value is 'linear'.
%
%     'hrtf',hrtf      HRTF database. This have to be in the TU-Berlin
%                      mat-format, see:
%                      https://dev.qu.tu-berlin.de/projects/measurements/wiki/IRs_file_format
%                      Default HRTF set is the 3m one from TU-Berlin measured
%                      with the KEMAR.
%
%     'lookup',lookup  Lookup table to map ITD values to angles. This can be
%                      created by the itd2anglelookuptable function. Default
%                      value is the lookup table
%                      wierstorf2013itd2anglelookup.mat that comes with AMT.
%
%
%   For the simulation of the wave field synthesis or stereophony setup this
%   functions depends on the Sound-Field-Synthesis Toolbox, which is available
%   here: http://github.com/sfstoolbox/sfs. It runs under Matlab and Octave. The
%   revision used to generate the figures in the corresponding paper is
%   a8914700a4.
%
%   See also: wierstorf2013estimateazimuth, dietz2011
%
%   References:
%     M. Dietz, S. D. Ewert, and V. Hohmann. Auditory model based direction
%     estimation of concurrent speakers from binaural signals. Speech
%     Communication, 53(5):592-605, 2011. [1]http ]
%     
%     H. Wierstorf, M. Geier, A. Raake, and S. Spors. A free database of
%     head-related impulse response measurements in the horizontal plane with
%     multiple distances. In Proceedings of the 130th Convention of the Audio
%     Engineering Society, 2011.
%     
%     H. Wierstorf, A. Raake, and S. Spors. Binaural assessment of
%     multi-channel reproduction. In J. Blauert, editor, The technology of
%     binaural listening, chapter 10. Springer, Berlin-Heidelberg-New York
%     NY, 2013.
%     
%     References
%     
%     1. http://www.sciencedirect.com/science/article/pii/S016763931000097X
%     
%
%   Url: http://amtoolbox.sourceforge.net/doc/binaural/wierstorf2013.php

% Copyright (C) 2009-2014 Peter L. Søndergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% AUTHOR: Hagen Wierstorf

%   Copyright (c) 2013   Assessment of IP-based Applications
%                        Technische Universitaet Berlin
%                        Ernst-Reuter-Platz 7, 10587 Berlin, Germany


%% ===== Checking of input parameters and dependencies ===================
nargmin = 7;
nargmax = 17;
error(nargchk(nargmin,nargmax,nargin));

definput.flags.method = {'stereo','wfs'};
definput.keyvals.array = 'linear';
definput.keyvals.nls = 2;
definput.keyvals.resolution = 21;
definput.keyvals.hrtf = [];
definput.keyvals.lookup = [];
definput.keyvals.showprogress = 1;
[flags,kv] = ...
    ltfatarghelper({'resolution','nls','array','hrtf','lookup','showprogress'},definput,varargin);
array = kv.array;
resolution = kv.resolution;
number_of_speakers = kv.nls;
hrtf = kv.hrtf;
lookup = kv.lookup;
showprogress = kv.showprogress;

% Checking for the Sound-Field-Synthesis Toolbox
if ~which('SFS_start')
    error(['%s: you need to install the Sound-Field-Synthesis Toolbox.\n', ...
        'You can download it at https://github.com/sfstoolbox/sfs.\n', ...
        'You need version 0.2.4 of the Toolbox (commit afe5c14359).'], ...
        upper(mfilename));
end

% Checking if we have only one position or if we have a whole listening area
if length(X)==1 && length(Y)==1
    resolution = 1;
    X = [X X];
    Y = [Y Y];
elseif length(X)==1
    X = [X X];
elseif length(Y)==1
    Y = [Y Y];
end


%% ===== Configuration ===================================================
% The following settings are all for the Sound Field Synthesis-Toolbox
% Binaural settings
% length of impulse response; this has two influences: 
% 1) longer impulse responses lead to a longer running time of the model
% 2) shorter impulse responses will not work if your WFS setup needs really long
% time shifting of single driving signals or you use HRTF with room reflections
% that rest longer than conf.N samples.
conf.N = 1024;
% Use no headphone compensation because we are not trying to listening to the
% signal
conf.ir.usehcomp = false;
conf.ir.hcompfile = '';
conf.ir.useinterpolation = true;
conf.ir.speechfile = '';
conf.ir.cellofile = '';
conf.ir.castanetsfile = '';
conf.ir.noisefile = '';
conf.ir.pinknoisefile = '';
%
% WFS settings
% dimensionality of the setup
conf.dimension = '2.5D';
% driving functions
conf.driving_functions = 'default';
% Use a WFS pre-equalization filter and specify its start frequency.
% Note, that the stop frequency will be calculated later with the aliasing
% frequency of your WFS setup.
conf.wfs.usehpre = true;
conf.wfs.hpretype = 'FIR';
conf.wfs.hpreflow = 50;
% Tapering window
conf.usetapwin = true;
conf.tapwinlen = 0.3;
% misc settings
conf.debug = 0;
conf.c = 343;
conf.fs = 44100;
conf.usefracdelay = 0;
conf.fracdelay_method = '';


%% ===== Loading of additional data ======================================
% Load default 3m TU-Berlin KEMAR HRTF from the net if no one is given to the
% function
if isempty(hrtf)
    % load HRTFs, see:
    % https://dev.qu.tu-berlin.de/projects/measurements/wiki/2010-11-kemar-anechoic
    [~,path] = download_hrtf('wierstorf2011_3m');
    load([path 'wierstorf2011_3m.mat']);
    hrtf = irs;
end
% Get sampling rate from the HRTFs
fs = hrtf.fs;
% Load lookup table from the AMT if no one is given to the function
if isempty(lookup)
    % load lookup table to map ITD values of the model to azimuth angles.
    % the lookup table was created using the same HRTF database
    path = which('amtstart');
    lookup = load([path(1:end-10) 'modelstages/wierstorf2013itd2anglelookup.mat']);
end


%% ===== Simulate the binaural ear signals ===============================
% Simulate a stereo setup
conf.secondary_sources.geometry = array;
% center of array
conf.secondary_sources.center = [0 0 0];
% initialize empty array
conf.secondary_sources.x0 = [];
% number of loudspeakers
conf.secondary_sources.number = number_of_speakers;
% length of array
conf.secondary_sources.size = L;
% get loudspeaker positions
x0 = secondary_source_positions(conf);
% selection of loudspeakers for WFS
if flags.do_wfs && strcmpi('circle',array)
    x0 = secondary_source_selection(x0,xs,src);
end

% calculate the stop frequency for the WFS pre-equalization filter
conf.wfs.hprefhigh = aliasing_frequency(x0,conf);

% get a grid of the listening positions
conf.resolution = resolution;
[~,~,x,y] = xy_grid(X,Y,conf);
% simulate the binaural impulse response
perceived_direction = zeros(length(x),length(y));
desired_direction = zeros(length(x),length(y));
localization_error = zeros(length(x),length(y));
for ii=1:length(x)
    if showprogress progressbar(ii,length(x)), end
    for jj=1:length(y)
        X = [x(ii) y(jj) 0];
        conf.xref = X;
        if flags.do_stereo
            % first loudspeaker
            ir1 = ir_point_source(X,phi,x0(1,1:3),hrtf,conf);
            % second loudspeaker
            ir2 = ir_point_source(X,phi,x0(2,1:3),hrtf,conf);
            % sum of both loudspeakers
            ir = (ir1+ir2)/2;
        else % WFS
            ir = ir_wfs(X,phi,xs,src,hrtf,conf);
        end
        % generate a 0.1s noise signal
        sig_noise = noise(fs/10,1,'white');
        % convolve with impulse response
        sig = auralize_ir(ir,sig_noise,1,conf);


        %% ===== Estimate the direction of arrival for the listener ==============
        % this is done by calculating ITDs with the dietz2011 binaural model, which are
        % then mapped to azimuth values with a lookup table
        %
        % estimate the perceived direction of arrival
        perceived_direction(ii,jj) = wierstorf2013estimateazimuth(sig,lookup,'dietz2011',0);
        % calculate the desired direction
        desired_direction(ii,jj) = source_direction(X,phi,xs,src);
        % calculate the localization error as the difference of both
    end
end
localization_error = perceived_direction - desired_direction;

end % of main function


%% ----- Subfunctions ----------------------------------------------------
function direction = source_direction(X,phi,xs,src)
    if strcmp('pw',src)
        [direction,~,~] = cart2sph(xs(1),xs(2),0);
        direction = deg(direction+phi);
    elseif strcmp('ps',src)
        x = xs-X;
        [direction,~,~] = cart2sph(x(1),x(2),0);
        % FIXME: this is not working with all points at the moment
        % For example place a stereo source at (0,0) and the listener at (0,-2)
        direction = deg(direction-phi);
    end
end

