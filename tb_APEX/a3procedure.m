function result=a3procedure(parameters, trials,type, id)
% function result=a3procedure(parameters, trials,type, id)
%
% 1. Description:
%       type 1 = adaptive procedure 
%       type 2 = constant procedure
%       type 3 = training procedure
%
% 2. Stand-alone example:
%       params.presentations = 1; % <presentations>1</presentations>
%       params.order = 'sequential'; % <order>sequential</order>
%       params.choices = 1; %<choices>1</choices>
%       params.nUp = 1; % <nUp>1</nUp>
%       params.nDown = 1; % <nDown>1</nDown>
%       params.adapt_parameter = 'itd'; % <adapt_parameter>itd</adapt_parameter>
%       params.start_value = 8; % <start_value>8</start_value>
%       params.stop_after_type = 'reversals'; % <stop_after_type>reversals</stop_after_type>
%       params.stop_after = 2; % <stop_after>2</stop_after>
%       params.min_value = -8; % <min_value>-8</min_value>
%       params.max_value = 8;  % <max_value>8</max_value>
%       params.larger_is_easier = 'false'; % <larger_is_easier>false</larger_is_easier>
%       params.repeat_first_until_correct = 'true'; % <repeat_first_until_correct>true</repeat_first_until_correct>
%       params.stepsizes = sprintf('\n\t<stepsize begin="0" size="4"/>\n\t<stepsize begin="1" size="2"/>\n\t<stepsize begin="2" size="1"/>\n');
%       a3procedure(params,2,1,'procedure1')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<4)
   id='';
else
   id=[' id= "' id '" '];
end

if (isstruct(parameters))
    parameters=a3toxml(parameters);
end

lf=sprintf('\n');
tb=sprintf('\t');
if (type==1)
    result=['<procedure xsi:type="apex:adaptiveProcedureType"' id '>' lf ];
elseif (type==2)
    result=['<procedure xsi:type="apex:constantProcedureType"' id '>' lf ];
elseif (type==3)
    result=['<procedure xsi:type="apex:trainingProcedureType"' id '>' lf ];
end
result=[result '<parameters>' lf tb parameters lf '</parameters>' lf trials lf '</procedure>' lf];
