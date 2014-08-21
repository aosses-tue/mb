function newstruct = sortStruct(s,field)
    for i = 1:numel(s)
        cdata{i} = s(i).(field);
    end
    [B,IDX] = sort(cdata);
    newstruct = s(IDX);
end