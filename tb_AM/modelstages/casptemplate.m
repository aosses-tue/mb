function [template_nonorm,ir_reference,ir_target, template]=casptemplate(target,reference,modelname,modelpars)
% function [template_nonorm,ir_reference,ir_target, template]=casptemplate(target,reference,modelname,modelpars)
%
% 1. Description:
%       CASPTEMPLATE  Generate a template for the optimal detector
%
%  CASPTEMPLATE(target,reference,modelname,modelpars) generates the template
%  needed for the optimal detector. CASPTEMPLATE will run the model
%  specified by modelname on the signals stored in target and reference
%  and generate the template from this.
%
%  If target or reference is a matrix, each column will be considered a
%  signal, and averaging will be done. This is usefull for stochastic
%  signals.
%
%   Url: http://amtoolbox.sourceforge.net/doc/modelstages/casptemplate.php
%
% 2. Stand-alone example:
%   fs          = 44100; % sampling frequency of the waveforms insig1 and insig2supra
%   target      = insig2supra + insig1;
%   reference   = insig1;
%   [template,ir_reference] = casptemplate(target,reference,'dau1996preproc',{fs});
% 
% Copyright (C) 2009-2014 Peter L. Soendergaard and Piotr Majdak.
% This file is part of AMToolbox version 0.9.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
    modelpars={};
end;

ntargets    = size(   target,2);
nreferences = size(reference,2);

%% ----- Compute average internal representation of the targets
ir_target=feval(modelname,target(:,1),modelpars{:});

for ii=2:ntargets
  ir_target = ir_target + feval(modelname,target(:,ii),modelpars{:});
end;

ir_target=ir_target/ntargets;

%% ----- Compute average internal representation of the references
ir_reference=feval(modelname,reference(:,1),modelpars{:});

for ii=2:nreferences
  ir_reference = ir_reference + feval(modelname,reference(:,ii),modelpars{:});
end;

ir_reference=ir_reference/nreferences;

% Compute the template as the difference between the average representation 
% of the targets and references.
template_nonorm = ir_target - ir_reference;

%% ----- Normalise to compensate for the increase in level ----

% Normalise across all dimenstions of the internal representation.
if nargout >= 4
    template=template_nonorm/rms(template_nonorm(:));
end

%OLDFORMAT

end
