function r20150630_update_experiments(fileres)
% function r20150630_update_experiments(fileres)
%
% 1. Description:
%       Files needed:
%           bProcessDataFig3- Daniel1997_Fig3_multiprocedure01-AO.apr.xml
%           bProcessDataFig5- Daniel1997_Fig5-AO-1.apr.xml
%   
% 2. Stand-alone example:
%           % Linux example, Alejandro's computer:
%           resultfile = '/home/alejandro/Documenten/Documenten-TUe/02-Experiments/2015-APEX-Rodrigo/APEX_shared/experiment/fluctuation_strength_results/AM_tones-fm-results-AO.apr';
%           r20150630_update_experiments(resultfile);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 30/06/2015
% Last update on: 01/07/2015 % Update this date manually
% Last use on   : 01/07/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bProcessDataFig5    = 1;

if nargin == 0
    if ~isunix % then is Windows
        dir_results = 'D:\Documenten-TUe\02-Experiments\2015-APEX-Rodrigo\APEX_shared\';
    else
        dir_results = '/home/alejandro/Documenten/Documenten-TUe/02-Experiments/2015-APEX-Rodrigo/APEX_shared/';
    end
    fileres = [dir_results 'experiment' delim 'fluctuation_strength_results' delim 'AM_tones-fm-results-AO.apr'];
end

trials1 = [11:19]; % put here the trial numbers corresponding to procedure1
trials2 = [21:29]; % put here the trial numbers corresponding to procedure2
answers_options = 0:25:200;

fileres = {'/home/alejandro/Documenten/Documenten-TUe/02-Experiments/2015-APEX-Rodrigo/APEX_shared/experiment/fluctuation_strength_results/AM_tones-fm-results-AO-test.apr', ...
           '/home/alejandro/Documenten/Documenten-TUe/02-Experiments/2015-APEX-Rodrigo/APEX_shared/experiment/fluctuation_strength_results/AM_tones-fm-slider-results-AO-test.apr'};

results1 = r20150630_plot_FS_results(fileres{1});
results2 = r20150630_plot_FS_results(fileres{2});

grand1 = [results1.scores1; results2.scores1];
grand2 = [results1.scores2; results2.scores2];

Mean2 = mean(grand2);
factor = 100/Mean2(6);

figure;
errorbar(1:9, mean(grand1),std(grand1)), hold on
errorbar(1:9, mean(grand2)*factor,std(grand2),'r');

error('Continue here')

if bProcessDataFig5
    
    script = 'apexresult.xsl';
        
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
                error('')
                % s1.score(1,i1) = max(answers_options)-answers_options( idx );
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
                error('')
                % s2.score(1,i2) = max(answers_options)-answers_options( idx );
            end
            
            i2 = i2+1;
            
        end
        
        
    end
    n = length(s1.trial)/length(trials1);
    test_fmod = length(trials1);

    scores1 = zeros(n,test_fmod);

    for i = 1:length(trials1)

        idx = find(trials1(i)==s1.trial);
        scores1(:,i) = transpose( s1.score(idx) );

    end

    n = length(s2.trial)/length(trials2);
    test_fmod = length(trials2);

    scores2 = zeros(n,test_fmod);

    for i = 1:length(trials1)

        idx = find(trials2(i)==s2.trial);
        scores2(:,i) = transpose( s2.score(idx) );

    end
    
    test_fmod = [0 0.25 0.5 1 2 4 8 16 32];
    MeansStd1 = mean(scores1);
    DevsStd1  = std(scores1);
    % tmp11 = percentile(scores1,25);
    % tmp12 = percentile(scores1,75);
    % tmp1 = [tmp11; tmp12];
    MeansStd2 = mean(scores2);
    DevsStd2  = std(scores2);
    factor2normalise = 100/MeansStd2(6); % corresponding to ref fmod of 4 Hz
    
    figure;
    subplot(3,1,1)
    errorbar(1:length(test_fmod),MeansStd1,DevsStd1), hold on
    ylabel('Relative FS [%]')
    xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fmod)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fmod)
    title('using Standard1')
    
    subplot(3,1,2)
    errorbar(1:length(test_fmod),MeansStd2,DevsStd2), hold on
    ylabel('Relative FS [%]')
    xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fmod)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fmod)
    title('using Standard2')
    fmod_points = 1:length(test_fmod);
    
    subplot(3,1,3)
    errorbar(fmod_points-0.1,MeansStd1                 ,DevsStd1), hold on
    errorbar(fmod_points+0.1,MeansStd2*factor2normalise,DevsStd2,'r'), hold on
    ylabel('Relative FS [%]')
    xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fmod)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fmod)
    legend('Std1','Std2')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
