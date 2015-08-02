function results = r20150802_plot_FS_results_v2(fileres)
% function results = r20150802_plot_FS_results_v2(fileres)
%
% 1. Description:
%           Plots Fluctuation strength results obtained using APEX.
%           The experiments were conducted using a magnitude estimation 
%           paradigm using two references (standard1 and standard2) related
%           to procedure1 and procedure2 of the APEX experiment file. The
%           test conditions are 9 in total
% 
% 2. Stand-alone example:
%           % Linux example, Alejandro's computer:
%           resultfile = '/home/alejandro/Documenten/Documenten-TUe/09-Training+activities/Master-thesis/01-2015-Rodrigo/20150729-pilot-APEX/AM-fm-Rodrigo.apr';
%           r20150802_plot_FS_results_v2(resultfile);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Original file name: r20150630_plot_FS_results
% Created on    : 02/08/2015
% Last update on: 02/08/2015 % Update this date manually
% Last use on   : 02/08/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    % Default input in Alejandro's computers
    if ~isunix % then is Windows
        dir_results = 'D:\Documenten-TUe\02-Experiments\2015-APEX-Rodrigo\APEX_shared\';
    else
        dir_results = '/home/alejandro/Documenten/Documenten-TUe/02-Experiments/2015-APEX-Rodrigo/APEX_shared/';
    end
    fileres = [dir_results 'experiment' delim 'fluctuation_strength_results' delim 'AM_tones-fm-results-AO-test.apr'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: customise parameters to your own Magnitude-estimation experiment
trials1 = [11:19 110:111]; % put here the trial numbers corresponding to procedure1
trials2 = [21:29 210:211]; % put here the trial numbers corresponding to procedure2
answers_options = 0:25:200;
test_fmod = [0 0.25 0.5 1 2 4 8 16 32 64 128];
text_XLabel = 'f_m_o_d [Hz]';
text_YLabel = 'Relative FS [%]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Reading the APEX result file

% % Uncomment the following lines in case not result section can be obtained from APEX result file:
% script = 'apexresult.xsl';
% forcetransform = 0; % set to 1 in case results were not processed yet...
% results = a3getresults(fileres,script,forcetransform);

results = a3getresults(fileres);
s = a3parseresults(results);
    % s-fields:
    %       procedure;  trial;      stimulus;   correctanswer;  corrector;  useranswer
    %       procedure1; trial_12;   stimtest12; 1;              false;      6
    %       YES         YES         YES         YES              NO          YES
    
N = length(s)/4/2;
trials1 = trials1(1:N);
trials2 = trials2(1:N);
test_fmod = test_fmod(1:N);

i1 = 1;
i2 = 1;

for i = 1:length(s)

    tmp = ( regexp(s(1,i).trial,['\d+\.?\d*'],'match') );
    tmp = str2num( tmp{1} ); 

    if      strcmp(s(1,i).procedure,'procedure1')

        s1.useranswer(1,i1) = str2num( s(1,i).useranswer );
        s1.trial(1,i1)      = tmp;
        try
            s1.referencefirst(1,i1) = str2num( s(1,i).correctanswer )-1;
            idx = s1.useranswer(1,i1);
            if s1.referencefirst(1,i1) == 1
                s1.score(1,i1) = answers_options( idx );
            else
                error('')
            end
        catch
            s1.score(1,i1) = str2num( s(1,i).useranswer );
        end
        
        % if s1.referencefirst(1,i1) == 1
        %     s1.score(1,i1) = answers_options( idx );
        % elseif s1.referencefirst(1,i1) == -2
        %     error('')
        %     s1.score(1,i1) = s(1,i).useranswer;
        % end

        i1 = i1+1;

    end

    if      strcmp(s(1,i).procedure,'procedure2')

        s2.useranswer(1,i2) = str2num( s(1,i).useranswer );
        s2.trial(1,i2)      = tmp;
        try
            s2.referencefirst(1,i2) = str2num( s(1,i).correctanswer )-1;
            idx = s2.useranswer(1,i2);
            if s2.referencefirst(1,i2) == 1
                s2.score(1,i2) = answers_options( idx );
            else
                error('')
            end
        catch
            s2.score(1,i2) = str2num( s(1,i).useranswer );
        end
        i2 = i2+1;

    end

end
n = length(s1.trial)/length(trials1);
scores1 = zeros(n,length(trials1));

for i = 1:length(trials1)

    idx = find(trials1(i)==s1.trial);
    scores1(:,i) = transpose( s1.score(idx) );

end

n = length(s2.trial)/length(trials2);
scores2 = zeros(n,length(trials2));

for i = 1:length(trials1)

    idx = find(trials2(i)==s2.trial);
    scores2(:,i) = transpose( s2.score(idx) );

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3: Processing the data
MeansStd1 = mean(scores1);
DevsStd1  = std(scores1);
MeansStd2 = mean(scores2);
DevsStd2  = std(scores2);
factor2normalise = 100/MeansStd2(6); % corresponding to ref fmod of 4 Hz
scores1norm = scores1*1;
scores2norm = scores2*factor2normalise;
MeansCombined = mean([MeansStd1; MeansStd2*factor2normalise]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4: Plotting the data

if nargout == 0
    figure;
    subplot(2,1,1)
    errorbar(1:length(test_fmod),MeansStd1,DevsStd1), hold on
    ylabel(text_YLabel) % ylabel('Relative FS [%]')
    xlabel(text_XLabel) % xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fmod)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fmod)
    title('FS using Standard1')

    subplot(2,1,2)
    errorbar(1:length(test_fmod),MeansStd2,DevsStd2), hold on
    ylabel(text_YLabel) % ylabel('Relative FS [%]')
    xlabel(text_XLabel) % xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fmod)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fmod)
    title('FS using Standard2')
    fmod_points = 1:length(test_fmod);

    figure;
    errorbar(fmod_points-0.1,MeansStd1                 ,DevsStd1), hold on
    errorbar(fmod_points+0.1,MeansStd2*factor2normalise,DevsStd2,'r'),
    plot(fmod_points,MeansCombined,'ko--','LineWidth',2)
    ylabel(text_YLabel) % ylabel('Relative FS [%]')
    xlabel(text_XLabel) % xlabel('fmod [Hz]')
    grid on
    xlim([0.5 length(test_fmod)+0.5])
    ha = gca;
    set(ha,'XTickLabel',test_fmod)
    legend('Std1','Std2','Combined')
end

results.raw = results;
results.scores1 = scores1;
results.scores2 = scores2;
results.scores1norm = scores1norm;
results.scores2norm = scores2norm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
