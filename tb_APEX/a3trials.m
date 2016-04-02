function result=a3trials(id_names, answers, id_screen, stimuli)
% function result=a3trials(names, answers, screens, stimuli)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (length(answers)==1)
        answers=repmat(answers(1), size(id_names,1),size(id_names,2));
end
if (~iscell(id_screen))
        id_screen=repmat(id_screen, size(id_names,1),size(id_names,2));
        error('FIXME: Convert to cell');
end

if (length(id_names)~=length(answers) && length(answers)~=0)
    error('Length of answers incorrect');
end
if (length(id_names)~=length(id_screen) )
    error('Length of screens incorrect');
end
if (length(id_names)~=length(stimuli))
    error('Length of stimuli incorrect');
end


lf=sprintf('\n');
tb=sprintf('\t');
result='<trials>';

for iF=1:length(id_names)
         
   temp=a3trial(id_names{iF}, id_screen{iF}, stimuli{iF}, answers{iF});
   
   result=[result temp];
end

result=[result '</trials>'];

return;

