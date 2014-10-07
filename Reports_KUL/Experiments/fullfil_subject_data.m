function pLT = fullfil_subject_data(pLT, subject)

nSubject = length(subject);
idxS1 = [];

for i = 1:nSubject
    idxS1 = [idxS1 length(find(pLT.aceData(:,pLT.column_subject)==subject(i)))];
end

idx_max = max(idxS1);

for j = 1:nSubject
    for i = 1:1:idx_max-idxS1(j)
        pLT.aceData = [pLT.aceData; [nan(1,size(pLT.aceData,2)-1) subject(j)]];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxS1 = [];

for i = 1:nSubject
    idxS1 = [idxS1 length(find(pLT.f0mData(:,pLT.column_subject)==subject(i)))];
end

idx_max = max(idxS1);

for j = 1:nSubject
    for i = 1:1:idx_max-idxS1(j)
        pLT.f0mData = [pLT.f0mData; [nan(1,size(pLT.f0mData,2)-1) subject(j)]];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idxS1 = [];

for i = 1:nSubject
    idxS1 = [idxS1 length(find(pLT.aceDataOwn(:,pLT.column_subject)==subject(i)))];
end

idx_max = max(idxS1);

for j = 1:nSubject
    for i = 1:1:idx_max-idxS1(j)
        pLT.aceDataOwn = [pLT.aceDataOwn; [nan(1,size(pLT.aceDataOwn,2)-1) subject(j)]];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end