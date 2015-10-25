function [succes,basepath] = download_hrtf(database)
%DOWNLOAD_HRTF download freely available HRTF data sets and store them in hrtf/data/
%
%
%   Url: http://amtoolbox.sourceforge.net/doc/hrtf/download_hrtf.php

% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%
% download_hrtf('wierstorf2011_3m');
% download_hrtf('wierstorf2011_2m');
% download_hrtf('wierstorf2011_1m');
% download_hrtf('wierstorf2011_0_5m');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ===== Database urls ===================================================
database_url = { ...
    'wierstorf2011_3m',   'https://dev.qu.tu-berlin.de/projects/measurements/repository/raw/2010-11-kemar-anechoic/mat/QU_KEMAR_anechoic_3m.mat';...
    'wierstorf2011_2m',   'https://dev.qu.tu-berlin.de/projects/measurements/repository/raw/2010-11-kemar-anechoic/mat/QU_KEMAR_anechoic_2m.mat';...
    'wierstorf2011_1m',   'https://dev.qu.tu-berlin.de/projects/measurements/repository/raw/2010-11-kemar-anechoic/mat/QU_KEMAR_anechoic_1m.mat';...
    'wierstorf2011_0_5m', 'https://dev.qu.tu-berlin.de/projects/measurements/repository/raw/2010-11-kemar-anechoic/mat/QU_KEMAR_anechoic_0.5m.mat';...
};


%% ===== Checking of input parameters ====================================


%% ===== Downloading the databases =======================================
% get the path to hrtf/
hpath = which('hrtfinit');  % find local path of hrtf repository
hpath = hpath(1:end-10);
basepath = [hpath filesep 'wierstorf2013'];
% create directory if it is not existing
if ~exist(basepath,'dir')
    mkdir(basepath);
end
for ii=1:size(database_url,1)
    if strcmpi(database_url{ii,1},database)
        % check if file is already downloaded
        filename = [basepath database_url{ii,1} '.mat'];
        if ~exist(filename,'file')
            succes = urlwrite(database_url{ii,2},filename);
        else
            succes = 1;
        end
        % end the loop if desired data base was stored
        return;
    end
end

end