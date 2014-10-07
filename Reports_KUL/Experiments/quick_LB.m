function [LB_value, AvgLB] = quick_LB(filename1, filename2)
% function [LB_value, AvgLB] = quick_LB(filename1, filename2)
%
% Example 1:
%       quick_LB; % then you have to select an LB APEX experiment file
%
% Example 2:
%       [LB_value, AvgLB] = quick_LB('LB_Com_104_UW-m5-ACE-AO.apr','LB_Com_104_UW-m5-F0m-AO.apr')
%       [LB_value, AvgLB] = quick_LB('LB_Com_294_UW-m20-ACE-AO.apr','LB_Com_294_UW-p0-F0m-AO.apr')
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    [filename_part1, filename_part2] = uigetfile('*',' LB Select an APEX experiment result file');
    filename1 = [filename_part2, filename_part1];
    
    bContinueLoading = input('Press 1 to continue loading LB Experiments: ');
    
    if bContinueLoading == 1
        [filename_part1, filename_part2] = uigetfile('*',' LB Select an APEX experiment result file');
        filename2 = [filename_part2, filename_part1];
    end
end 

[LB_value, ref_value] = get_balance_levels(filename1);

if nargin == 2
    [LB_value(2,:), ref_value(2,:)] = get_balance_levels(filename2);
end

AvgLB = mean(LB_value);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end