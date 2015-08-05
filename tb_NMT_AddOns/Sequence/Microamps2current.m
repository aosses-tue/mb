function CU = Microamps2current(uA,ci_type)
% function CU = Microamps2current(uA, ci_type)

if ~exist('ci_type','var')
    ci_type = 'CIC3';
end

if strcmp(ci_type,'CIC3')
    CU = 255/log(175)*log(uA/10);
end

