function s = varargin_to_struct(varargin)

s = struct;

properties = varargin;
while length(properties) >= 2,
   prop = properties{1};
   val = properties{2};
   s.(prop) = val;
   properties = properties(3:end);
end

if length(properties) == 1
    n = fieldnames(properties{1});
    for k = 1:length(n)
        if ~isfield(s, n{k})
            s.(n{k}) = properties{1}.(n{k});
        end
    end
end
