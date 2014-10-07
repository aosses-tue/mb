function [LB_value, ref_value] = get_balance_levels(base_filename, InitialsSubject)
% [LB_value, ref_value] = function get_balance_levels(base_filename, InitialsSubject)
%
% This function get the results of the Balance Loudness experiments for the 
% specified base_filename (APEX result experiment) which matches with InitialSubject
%
% colResult = 6 (automate this in the future): column where the results are (APEX format)
% numLastSamples 6 (automate this in the future)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
    Final_name = base_filename;
else
    Final_name = dir([cd delim base_filename '*' InitialsSubject '*.apr']);
    if length(Final_name)==0
        LB_value = NaN;
        ref_value = 0;
        return;
    end
    Final_name = Final_name.name;
end
        
[resultsm, parametersm, generalm] = a3getresults(Final_name);
        
colResult       = 6;
Result          = 0;
numLastSamples  = 6;

k               = 0;

for i = 1:numLastSamples
    tempResult = strsplit(resultsm{end-1-i},';');
    try
        Result = Result + str2num( tempResult{colResult} ); 
    catch
        Result = Result + str2num( tempResult{3} ); 
    end
end

LB_value = Result / numLastSamples;
ref_value = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end