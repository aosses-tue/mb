function s = ef(s, fieldName, defaultVal)
% function s = ef(s, fieldName, defaultVal)
%
% 1. Description:
%       Same than Ensure_field, but without notifying when a field is assigned
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(s, fieldName)
    s = setfield(s, fieldName, defaultVal);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end