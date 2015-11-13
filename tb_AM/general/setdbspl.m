function inoutsig = setdbspl(inoutsig,lvl,varargin)
% function inoutsig = setdbspl(inoutsig,lvl,varargin)
% 
%	Set level of signal in dB
%   Convention (RMS):     100 dB SPL =   0 dBFS
%                          75        = -25 dBFS
%                          65        = -35 dBFS
% 
%   Usage: outsig = setdbspl(insig,lvl);
%          outsig = setdbspl(insig,lvl,'ac');
%
%   SETDBSPL(insig,lvl) sets the SPL (sound pressure level) of the signal
%   insig to lvl dB, using the convention that a pure tone with an RMS value
%   of 1 corresponds to 100 dB SPL.
%
%   SETDBSPL(lvl) returns a scaling constant that will scale a signal
%   with RMS value of 1 to the correct level. 
%
%   If the input is a matrix, it is assumed that each column is a signal.
%
%   SETDBSPL(insig,lvl,'ac') does the same, but considers only the AC
%   component of the signal (i.e. the mean is removed).
%
%   References:
%     B. C. J. Moore. An Introduction to the Psychology of Hearing. Academic
%     Press, 5th edition, 2003.
%
%   Url: http://amtoolbox.sourceforge.net/doc/general/setdbspl.php
%
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5-0.9.7
% 
%   See also: dbspl
% 
%   Author: Peter L. Soendergaard, 2009
%   Comments edited by: Alejandro Osses, a.osses@tue.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------ Checking of input parameters ---------

definput.flags.mean={'noac','ac'};
definput.keyvals.dboffset=100;
[flags,kv]=ltfatarghelper({'dboffset'},definput,varargin);

error(nargchk(1,5,nargin));

if ~isnumeric(inoutsig)
  error('%s: insig must be numeric.',upper(mfilename));
end;

% In the code below, "setdbspl" obtains the reference level from "dbspl"
% by calling "dbspl(1)", which will return only the offset measured in dB.

if nargin==1
  % Special mode, only the level has been given
  lvl=inoutsig;
  
  if ~isscalar(lvl) 
    error('%s: lvl must be a scalar.',upper(mfilename));
  end;

  inoutsig=gaindb(1,lvl-dbspl(1));
  return;
end;


if ~isnumeric(lvl) || ~isscalar(lvl) 
  error('%s: lvl must be a scalar.',upper(mfilename));
end;

% if (nargin<3) || (~ischar(options))
%   options='';
% end;


% ------ Computation --------------------------

if isvector(inoutsig)
  inoutsig = gaindb(inoutsig/rms(inoutsig,flags.mean),...
    lvl-dbspl(1,'dboffset',kv.dboffset));
else
  % If we have a matrix, set the level for every column.
  for ii=1:size(inoutsig,2);
    inoutsig(:,ii) = gaindb(inoutsig(:,ii)/rms(inoutsig(:,ii),flags.mean),...
      lvl-dbspl(1,'dboffset',kv.dboffset));
  end;
end;
