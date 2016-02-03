function results = Get_APEX_slider_results(fileres)
% function results = Get_APEX_slider_results(fileres)
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
%           resultfile = '/home/alejandro/Documenten/Documenten-TUe/02-Experiments/2015-APEX-Rodrigo/APEX_shared/experiment/fluctuation_strength_results/AM_tones-fm-results-AO-test.apr';
%           Get_APEX_slider_results(resultfile);
% 
%           resultfile = '/home/alejandro/Documenten/Documenten-TUe/09-Training+activities/Master-thesis/01-2015-Rodrigo/20150729-pilot-APEX/AM-fm-Rodrigo.apr';
%           Get_APEX_slider_results(resultfile);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses V., HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 30/06/2015
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
test_fmod = [0 0.25 0.5 1 2 4 8 16 32];
text_XLabel = 'f_m_o_d [Hz]';
text_YLabel = 'Relative FS [%]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Reading the APEX result file

% % Uncomment the following lines in case not result section can be obtained from APEX result file:
% script = 'apexresult.xsl';
% forcetransform = 0; % set to 1 in case results were not processed yet...
% results = a3getresults(fileres,script,forcetransform);

results = a3getresults(fileres);
s       = a3parseresults(results);
Nrep    = 4;

    % s-fields:
    %       procedure;  trial;      stimulus;   correctanswer;  corrector;  useranswer
    %       procedure1; trial_12;   stimtest12; 1;              false;      6
    %       YES         YES         YES         YES              NO          YES
i1 = 1;
i2 = 1;

for i = 1:length(s)

    tmp = ( regexp(s(1,i).trial,['\d+\.?\d*'],'match') );
    tmp = str2num( tmp{1} ); % gets trial number

    if strcmp(s(1,i).procedure,'procedure1')

        s1.useranswer(1,i1) = str2num( s(1,i).useranswer );
        s1.trial(1,i1)      = tmp;
        s1.score(1,i1)      = str2num( s(1,i).useranswer );
        trials1(i1)         = tmp;
        
        i1 = i1+1;

    end

    if strcmp(s(1,i).procedure,'procedure2')

        s2.useranswer(1,i2) = str2num( s(1,i).useranswer );
        s2.trial(1,i2)      = tmp;
        s2.score(1,i2)      = str2num( s(1,i).useranswer );
        trials2(i2)         = tmp;
        
        i2 = i2+1;

    end

end

trials1u = unique(sort(trials1));
trials2u = unique(sort(trials2));

scores1 = nan(Nrep,length(trials1u));
scores2 = nan(Nrep,length(trials2u));

for i = 1:length(trials1u)

    idx = find(trials1u(i)==s1.trial);
    Lidx = length(idx);
    scores1(1:Lidx,i) = transpose( s1.score(idx) );

end

for i = 1:length(trials2u)

    idx = find(trials2u(i)==s2.trial);
    Lidx = length(idx);
    scores2(1:Lidx,i) = transpose( s2.score(idx) );

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3: Processing the data
MeansStd1 = mean(scores1);
DevsStd1  = std(scores1);
MeansStd2 = mean(scores2);
DevsStd2  = std(scores2);
try
    factor2normalise = 100/MeansStd2(6); % corresponding to ref fmod of 4 Hz
    MeansCombined = mean([MeansStd1; MeansStd2*factor2normalise]);
catch
    warning('no normalisation factor stored');
end
    
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

results.scores1 = scores1;
results.scores2 = scores2;
try
    results.factor2normalise = factor2normalise;
end
results.raw = results;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
