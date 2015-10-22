function bp = amtbasepath;
% function bp = amtbasepath;
% 
% 1. Description:
% 
%AMTBASEPATH  The base path of the AMT installation
%   Usage: bp = amtbasepath;
%
%   AMTBASEPATH returns the top level directory in which the AMT
%   files are installed.
%
%   See also: amtstart
%
%   Url: http://amtoolbox.sourceforge.net/doc/amtbasepath.php
%
% Copyright (C) 2009-2015 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5-0.9.7
%
  
f=mfilename('fullpath');

bp = f(1:end-11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
