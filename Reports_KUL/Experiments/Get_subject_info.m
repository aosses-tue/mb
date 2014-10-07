function info = Get_subject_info(subject)
% function info = Get_subject_info(subject)
%
% % Example:
%       subject = 'ci-Maria_Brughmans';
%       info = Get_subject_info(subject);
% Programmed by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list_of_subjects = {'ci-Romain_Peeters'     ,11; ...
                    'ci-Maria_Brughmans'    ,12; ...
                    'ci-Patrick_Meul'       ,13; ...
                    'ci-Jan_leys'           ,14; ...
                    'ci-Wouter_David'       ,15; ...
                    'ci-Jean-Marie_Daumerie',16; ...
                    'ci-Julia_Schoolmeesters',17};

if strcmp( subject(end), delim )
    subject = subject(1:end-1);
end
                
for i = 1:size(list_of_subjects,1)
    if strcmp(subject,list_of_subjects{i,1})
        info.nSubject = list_of_subjects{i,2};
    end
end

tmpSubject = strsplit(subject,'-');
info.type       = tmpSubject{1};

tmpSubject = strsplit(subject,'_');
info.firstname  = tmpSubject{1};

info.Initials   = KeepCapitalLetters(subject); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end