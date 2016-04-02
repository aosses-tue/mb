function result=a3stimulus(stimid, datablocks,fixedparameters,variableparameters,simultaneous)
% function temp=a3stimulus(id, datablocks,fixedparameters,variableparameters,simultaneous)
%
%   1. Description:
%           a3stimulus(id, datablocks,fixedparameters,variableparameters,
%           simultaneous) make simple stimulus with n datablocks in parallel
%           datablocks is a cell of datablocks parameters is a struct of 
%           fixed parameters if simultaneous==1, datablocks are put in 
%           parallel, otherwise in series
% 
% 2. Stand-alone example:
%       
%       stimid = 'stim1';
%       datablocks = {'data_stim1','data_noise1'};
%       result = a3stimulus(stimid, datablocks); %, fixedparameters, variableparameters, simultaneous);
%       disp(result);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<5)
    simultaneous = 1;
end

if (nargin < 3)
    fixedparameters=struct;
end
if (nargin < 4)
    variableparameters=struct;
end

if (~iscell(datablocks))
    datablocks={datablocks};
end

result=['<stimulus id="' stimid '">' lf];
result=[result tb '<datablocks>' lf];
if size(datablocks,2)>1
    result=[result tb '<sequential>' lf];
    result=[result tb '<simultaneous>' lf];
else
    result=[result tb '<sequential>' lf];
end
for iF=1:length(datablocks)
    result=[result tb tb '<datablock id="' datablocks{iF} '"/>' lf];
end
if size(datablocks,2)>1
    result=[result tb '</simultaneous>' lf];
%     temp=[temp tb '<simultaneous>' lf];
%     for iF=1:size(datablocks,2)
%     temp=[temp tb tb '<datablock id="' datablocks{1,iF} '"/>' lf];
%     temp=[temp tb '</simultaneous>' lf];
%     end
    result=[result tb '</sequential>' lf];
else
    result=[result tb '</sequential>' lf];
end
result=[result tb '</datablocks>' lf];

result=[result tb '<variableParameters>' lf];
result=[result tb a3showparams(variableparameters) ];
result=[result tb '</variableParameters>' lf];

result=[result tb '<fixedParameters>' lf];
result=[result tb a3showparams(fixedparameters) ];
result=[result tb '</fixedParameters>' lf];

result=[result '</stimulus>' lf];
