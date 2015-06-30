function r20150630_update_experiments(dir_results)
% function r20150630_update_experiments(dir_results)
%
% 1. Description:
%       Files needed:
%           bProcessDataFig3- Daniel1997_Fig3_multiprocedure01-AO.apr.xml
%           bProcessDataFig5- Daniel1997_Fig5-AO-1.apr.xml
%   
%       The APEX result files should be inside the dir_results folder
% 
% 2. Stand-alone example:
%           % Windows example, Alejandro's computer:
%           dir_results = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Roughness\';
%           r20150306_update_experiments(dir_results);
% 
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/03/2015
% Last update on: 05/03/2015 % Update this date manually
% Last use on   : 05/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    if ~isunix % windows
        dir_results = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Roughness\';
    else % linux, mac
        dir_results = uigetdir(pwd,'Choose a directory where APEX result files could be found...');
        dir_results = [dir_results delim];
    end
end

bDiary = 0;
Diary(mfilename,bDiary);

bProcessDataFig3    = 1;
bProcessDataFig5    = 1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if bProcessDataFig3
%     
%     % fileres = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\Pilot1\Daniel1997_Fig3_multiprocedure01-AO.apr.xml';
%     % answers_options = [0 25 75 100 125 150 175 200];
%     % fileres = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\Daniel1997_Fig3_multiprocedure01-AO.apr.xml';
%     % answers_options = [0 25 50 75 100 125 150 175 200];
%     
%     fileres = [dir_results 'Daniel1997_Fig3_multiprocedure01-AO.apr.xml'];
%     answers_options = [0 25 50 75 100 125 150 175 200];
%     
%     script = 'apexresult.xsl';
%     
%     forcetransform = 1;
%     [results,parameters,general] = a3getresults(fileres,script,forcetransform)
%     s = a3parseresults(results);
%     % s-fields:
%     %       procedure;  trial;      stimulus;   correctanswer;  corrector;  useranswer
%     %       procedure1; trial_12;   stimtest12; 1;              false;      6
%     %       YES         YES         YES         YES              NO          YES
%     i1 = 1;
%     i2 = 1;
%     
%     for i = 1:length(s)
%         
%         tmp = ( regexp(s(1,i).trial,['\d+\.?\d*'],'match') );
%         tmp = str2num( tmp{1} ); 
%         
%         if      strcmp(s(1,i).procedure,'procedure1')
%             
%             s1.useranswer(1,i1) = str2num( s(1,i).useranswer );
%             s1.trial(1,i1)      = tmp;
%             s1.referencefirst(1,i1) = str2num( s(1,i).correctanswer )-1;
%             
%             idx = s1.useranswer(1,i1);
%             if s1.referencefirst(1,i1) == 1
%                 s1.score(1,i1) = answers_options( idx );
%             else
%                 s1.score(1,i1) = max(answers_options)-answers_options( idx );
%             end
%             
%             i1 = i1+1;
%             
%         elseif  strcmp(s(1,i).procedure,'procedure2')
%             
%             s2.useranswer(1,i2) = str2num( s(1,i).useranswer );
%             s2.trial(1,i2)      = tmp;
%             s2.referencefirst(1,i2) = str2num( s(1,i).correctanswer )-1;
%             
%             idx = s2.useranswer(1,i2);
%             if s2.referencefirst(1,i2) == 1
%                 s2.score(1,i2) = answers_options( idx );
%             else
%                 s2.score(1,i2) = max(answers_options)-answers_options( idx );
%             end
%             i2 = i2+1;
%         end
%     end
%     trials = [11:18];
%     n = length(s1.trial)/length(trials);
%     m = length(trials);
% 
%     scores1 = zeros(n,m);
%     scores2 = zeros(n,m);
% 
%     for i = 1:length(trials)
% 
%         idx = find(trials(i)==s1.trial);
%         scores1(:,i) = transpose( s1.score(idx) );
% 
%     end
% 
%     trials = [21:28];
%     for i = 1:length(trials)
% 
%         idx = find(trials(i)==s2.trial);
%         scores2(:,i) = transpose( s2.score(idx) );
% 
%     end
% 
%     fmod = [40:10:110];
%     close all
%     figure;
%     errorbar(fmod,mean(scores1),std(scores1)), hold on
%     title('Scores 1: ref at 1000 Hz')
%     xlabel('f_m_o_d [Hz]')
%     grid on
% 
%     figure;
%     errorbar(fmod,mean(scores2),std(scores2)), hold on
%     title('Scores 2: ref at 500 Hz')
%     xlabel('f_m_o_d [Hz]')
%     grid on
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileres = 'D:\Documenten-TUe\02-Experiments\2015-APEX-Rodrigo\APEX_shared\experiment\fluctuation_strength_results\AM_tones-fm-results-AO.apr';

if bProcessDataFig5
    
    % fileres = [dir_results 'Daniel1997_Fig5-AO-1.apr.xml'];
    script = 'apexresult.xsl';
    answers_options = 0:25:200;
    
    forcetransform = 1;
    [results,parameters,general] = a3getresults(fileres,script,forcetransform)
    s = a3parseresults(results);
    % s-fields:
    %       procedure;  trial;      stimulus;   correctanswer;  corrector;  useranswer
    %       procedure1; trial_12;   stimtest12; 1;              false;      6
    %       YES         YES         YES         YES              NO          YES
    i1 = 1;
    i2 = 1;
    
    for i = 1:length(s)
        
        tmp = ( regexp(s(1,i).trial,['\d+\.?\d*'],'match') );
        tmp = str2num( tmp{1} ); 
        
        if      strcmp(s(1,i).procedure,'procedure1')
            
            s1.useranswer(1,i1) = str2num( s(1,i).useranswer );
            s1.trial(1,i1)      = tmp;
            s1.referencefirst(1,i1) = str2num( s(1,i).correctanswer )-1;
            
            idx = s1.useranswer(1,i1);
            if s1.referencefirst(1,i1) == 1
                s1.score(1,i1) = answers_options( idx );
            else
                s1.score(1,i1) = max(answers_options)-answers_options( idx );
            end
            
            i1 = i1+1;
            
        end
        
        if      strcmp(s(1,i).procedure,'procedure2')
            
            s2.useranswer(1,i2) = str2num( s(1,i).useranswer );
            s2.trial(1,i2)      = tmp;
            s2.referencefirst(1,i2) = str2num( s(1,i).correctanswer )-1;
            
            idx = s2.useranswer(1,i2);
            if s2.referencefirst(1,i2) == 1
                s2.score(1,i2) = answers_options( idx );
            else
                s2.score(1,i2) = max(answers_options)-answers_options( idx );
            end
            
            i2 = i2+1;
            
        end
        
        
    end
    trials1 = [11:19];
    n = length(s1.trial)/length(trials1);
    m = length(trials1);

    scores1 = zeros(n,m);

    for i = 1:length(trials1)

        idx = find(trials1(i)==s1.trial);
        scores1(:,i) = transpose( s1.score(idx) );

    end

    trials2 = [21:29];
    n = length(s2.trial)/length(trials2);
    m = length(trials2);

    scores2 = zeros(n,m);

    for i = 1:length(trials1)

        idx = find(trials2(i)==s2.trial);
        scores2(:,i) = transpose( s2.score(idx) );

    end
    
    m = [0 0.25 0.5 1 2 4 8 16 32];
    % Rest = 100*( m.^1.6 );
    figure;
    subplot(2,1,1)
    errorbar(1:length(m),mean(scores1),std(scores1)), hold on
    % plot(100*m,Rest,'r--')
    ylabel('Relative FS [%]')
    xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(m)+0.5])
    ha = gca;
    set(ha,'XTickLabels',m)
    
    subplot(2,1,2)
    errorbar(1:length(m),mean(scores2),std(scores1)), hold on
    % plot(100*m,Rest,'r--')
    ylabel('Relative FS [%]')
    xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(m)+0.5])
    ha = gca;
    set(ha,'XTickLabels',m)
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
