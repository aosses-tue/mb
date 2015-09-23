function ir_reference = casprepresentation(referencestim,modelname,modelpars)
% function ir_reference = casprepresentation(referencestim,modelname,modelpars)
%
% 1. Description:
%       CASPREPRESENTATION  Generate an internal representation to be used
%       by the optimal detector
%
%  CASPTEMPLATE(target,reference,modelname,modelpars) generates the template
%  needed for the optimal detector. CASPTEMPLATE will run the model specified 
%  by modelname on the signals stored in target and reference and generate 
%  the template from this.
%
%  If target or reference is a matrix, each column will be considered a
%  signal, and averaging will be done. This is usefull for stochastic
%  signals.
%
%   Url: http://amtoolbox.sourceforge.net/doc/modelstages/casptemplate.php
%
% 2. Stand-alone example:
%   fs          = 44100; % sampling frequency of the waveforms insig1 and insig2supra
%   target      = insig1;
%   reference   = insig2supra;
%   [template,ir_reference] = casptemplate(target,reference,'dau1996preproc',{fs});
% 
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3
    modelpars={};
end;

nreferences = size(referencestim,2);

%% ----- Compute average internal representation of the references
ir_reference=feval(modelname,referencestim(:,1),modelpars{:});

for ii=2:nreferences
  ir_reference = ir_reference + feval(modelname,referencestim(:,ii),modelpars{:});
end;

ir_reference=ir_reference/nreferences;

%OLDFORMAT

