function r20150306_update_experiments
% function r20150306_update_experiments
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014
% Created on    : 05/03/2015
% Last update on: 05/03/2015 % Update this date manually
% Last use on   : 05/03/2015 % Update this date manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bDiary = 0;
Diary(mfilename,bDiary);

bDoRoughFig3        = 0;
bDoRoughFig5        = 0;
bDoFluct            = 1;
bProcessDataFig3    = 1;
bProcessDataFig5    = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bDoRoughFig3

    % Daniel, Fig.3(c)
    lvl = 60;
    fc  = [125 250 500 1000 2000 4000 8000];
    fmod = 40:10:120;
    m_or_d = 100;
    option = 'm';
    fs = 44100;
    dur = 0.8;
    dur_ramp_ms = 20;
    start_phase = 0;
    
    for j = 1:length(fc)
        for i = 1:length(fmod)
            tmp = Generate_test_AMtone(fc(j),fmod(i),m_or_d,option,lvl,fs,dur,dur_ramp_ms,start_phase);
        end
    end
    
end

if bDoRoughFig5

    % Daniel, Fig.5
    lvl = 70;
    fc  = 1000;
    fmod = 70;
    m_or_d = 0:4:100;
    option = 'm';
    fs = 44100;
    dur = 0.8;
    dur_ramp_ms = 20;
    start_phase = 0;
    
    for j = 1:length(m_or_d)
        tmp = Generate_test_AMtone(fc,fmod,m_or_d(j),option,lvl,fs,dur,dur_ramp_ms,start_phase);
    end
    
    % This is for an additional figure...
    fc  = 500;
    
    for j = 1:length(m_or_d)
        tmp = Generate_test_AMtone(fc,fmod,m_or_d(j),option,lvl,fs,dur,dur_ramp_ms,start_phase);
    end
end

if bDoFluct

    lvl = 60;
    fc  = 1000;
    fmod = [0.5 1 2 4 8 16 32 50 70 90 110 130 150];
    m_or_d = 100;
    option = 'm';
    fs = 44100;
    dur = 4;
    dur_ramp_ms = 20;
    start_phase = 0;
    
    for j = 1:length(fmod)
        tmp = Generate_test_AMtone(fc,fmod(j),m_or_d,option,lvl,fs,dur,dur_ramp_ms,start_phase);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bProcessDataFig3
    
    % fileres = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\Pilot1\Daniel1997_Fig3_multiprocedure01-AO.apr.xml';
    % answers_options = [0 25 75 100 125 150 175 200];
    % fileres = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\Daniel1997_Fig3_multiprocedure01-AO.apr.xml';
    % answers_options = [0 25 50 75 100 125 150 175 200];
    
    fileres = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\Daniel1997_Fig3_multiprocedure01-AO.apr.xml'
    answers_options = [0 25 50 75 100 125 150 175 200];
    
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
                s1.score(1,i1) = max(answers_options)-answers_options( idx );
            end
            
            i1 = i1+1;
            
        elseif  strcmp(s(1,i).procedure,'procedure2')
            
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
    trials = [11:18];
    n = length(s1.trial)/length(trials);
    m = length(trials);

    scores1 = zeros(n,m);
    scores2 = zeros(n,m);

    for i = 1:length(trials)

        idx = find(trials(i)==s1.trial);
        scores1(:,i) = transpose( s1.score(idx) );

    end

    trials = [21:28];
    for i = 1:length(trials)

        idx = find(trials(i)==s2.trial);
        scores2(:,i) = transpose( s2.score(idx) );

    end

    fmod = [40:10:110];
    close all
    figure;
    errorbar(fmod,mean(scores1),std(scores1)), hold on
    title('Scores 1: ref at 1000 Hz')
    xlabel('f_m_o_d [Hz]')
    grid on

    figure;
    errorbar(fmod,mean(scores2),std(scores2)), hold on
    title('Scores 2: ref at 500 Hz')
    xlabel('f_m_o_d [Hz]')
    grid on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bProcessDataFig5
    
    fileres = 'D:\Documenten-TUe\02-Experiments\2015-APEX-my-experiments\Experiments\Roughness\Daniel1997_Fig5-AO-1.apr.xml'
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
    end
    trials = [10:20];
    n = length(s1.trial)/length(trials);
    m = length(trials);

    scores1 = zeros(n,m);

    for i = 1:length(trials)

        idx = find(trials(i)==s1.trial);
        scores1(:,i) = transpose( s1.score(idx) );

    end

    m = [0:.1:1];
    % Rest = 100*( m.^1.6 );
    figure;
    errorbar(100*m,mean(scores1),std(scores1)), hold on
    % plot(100*m,Rest,'r--')
    ylabel('Relative R [%]')
    xlabel('modulation index m [%]')
    grid on
    xlim([0 100])
end

if bDiary
	diary off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['EOF: ' mfilename '.m'])
