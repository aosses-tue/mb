function [race, rf0m, stRove] = plot_PR_results( filename, InitialSubject, ref_note, subject, do_plot, bTruncate)
% function [race, rf0m, stRove] = plot_PR_results( filename, InitialSubject, ref_note, subject, do_plot, bTruncate)
%
% Assumption: 4 semitones being compared to the reference
%
% INPUTS:
% result_file: cell_array containing xml result files usually in order ACE, F0mod
%
% Outputs:
%   race
%   rf0m
%   stRove - structure with Roving information
%
% Programmed by Alejandro, adapted from a Matthias' script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('bTruncate','var')
    bTruncate = 0;
end

if nargin ~= 1
    if bTruncate == 0
        result_files = {[filename 'ACE-' InitialSubject '.apr'], ...
                        [filename 'F0m-' InitialSubject '.apr']};
    else
        result_files = {[filename 'ACE-' InitialSubject '-trunc.apr'], ...
                        [filename 'F0m-' InitialSubject '-trunc.apr']};
    end
    txtLegend =  {'ACE', 'F0mod'};
else
    result_files = {filename};
    txtLegend   = {'Instant result'};
    % subject     = ': Current subject';
    do_plot     = 1;
end

if nargin < 3
    ref_note = '';
end

stRove = [];

% TimesPresented: number of times in which each pair reference - tone was presented (per trial)
try
    [h, xx, xx, corrector, TimesPresented, roving] = eval_results(result_files, 0, txtLegend);
catch
    if bTruncate == 0
        result_files = {[filename(1:end-1) '_ACE-LB-ACE-' InitialSubject '.apr'], ... % Loudness balanced experiments
                        [filename(1:end-1) '_F0m-LB-F0m-' InitialSubject '.apr']};
    else
        result_files = {[filename(1:end-1) '_ACE-LB-ACE-' InitialSubject '-trunc.apr'], ... % Loudness balanced experiments
                        [filename(1:end-1) '_F0m-LB-F0m-' InitialSubject '-trunc.apr']};
    end
	[h, xx, xx, corrector, TimesPresented, roving] = eval_results(result_files, 0, txtLegend);
end

hace        = h{1}{1};

if nargout >= 2
    hf0m        = h{2}{1};
end

race        = zeros(1,4);
rf0m        = zeros(1,4);

switch ref_note
    
    case {'Gsh2','G#2'}
        xticks = {'A2', 'Ash2', 'B2', 'C3'};
    case 'C3'
        xticks = {'Csh3', 'D3', 'Dsh3', 'E3'};
    case 'E3'
        xticks = {'F3', 'Fsh3', 'G3', 'Gsh3'};
    case {'Gsh3','G#3'}
        xticks = {'A3', 'Ash3', 'B3', 'C4'};
    case 'C4'
        xticks = {'Csh4', 'D4', 'Dsh4', 'E4'};
    otherwise
        xticks = {'1','2','3','4'};
end        

xticks_in_semitones = {'1','2','3','4'};

if xticks{1} ~= '1'

    for i=1:length(xticks)
        
        if length(hace{1,1}) ~= 1
            ToCompare = xticks;
        else
            ToCompare = xticks_in_semitones;
        end
        
        for j=1:size(hace,1)
            if( strcmp(hace{j,1}, ToCompare{i} ))
                race(i) = hace{j, 2}/TimesPresented*100; % percentage
            end
        end

        for j=1:size(hf0m, 1)
            if( strcmp(hf0m{j,1}, ToCompare{i} ))
                rf0m(i) = hf0m{j, 2}/TimesPresented*100; % percentage
            end
        end

    end
    
else
    
    bFound = 0;
    count_possibilities = 1;
    
    while bFound == 0
        
        if length(hace{1,1}) > 1 % then histogram returned semitones in format 'C_2'
            ref_note_possible = {'Gsh2','C3','E3','Gsh3','C4'};
            ref_note = ref_note_possible{count_possibilities};

            switch ref_note

                case {'Gsh2','G#2'}
                    xticks = {'A2', 'Ash2', 'B2', 'C3'};
                case 'C3'
                    xticks = {'Csh3', 'D3', 'Dsh3', 'E3'};
                case 'E3'
                    xticks = {'F3', 'Fsh3', 'G3', 'Gsh3'};
                case {'Gsh3','G#3'}
                    xticks = {'A3', 'Ash3', 'B3', 'C4'};
                case 'C4'
                    xticks = {'Csh4', 'D4', 'Dsh4', 'E4'};

            end   
        end
        
        for i=1:length(xticks)

            for j=1:size(hace, 1)
                if( strcmp(hace{j,1}, xticks{i} ))
                    bFound = 1;
                    race(i) = hace{j, 2}/TimesPresented*100; % percentage
                end
            end

            if nargout >= 2
                for j=1:size(hf0m, 1)
                    if( strcmp(hf0m{j,1}, xticks{i} ))
                        rf0m(i) = hf0m{j, 2}/TimesPresented*100; % percentage
                    end
                end
            end

        end
        
        count_possibilities = count_possibilities + 1;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roving analysis and plotting
if length(roving) ~= 0
    % Roving limits:
    try
        minmax_rove = minmax(roving{1}');
        min_rove = minmax_rove(1);
        max_rove = minmax_rove(2);

        vls_rove = [-5:1:5];
        idx_rove = 1:length(vls_rove);
        stRove.roveACEScore(idx_rove) = NaN;
        stRove.roveF0mScore(idx_rove) = NaN;
        dist2minrove = min_rove - min(vls_rove) + 1;

        for i = 1:length(roving)
            for j = round(min_rove):round(max_rove)-1
                idx = find( floor(roving{i})==j );
                Avg_ans(i,j-round(min_rove)+1) = mean(corrector{i}(idx));
            end
            
            [Avg_num x_dB] = hist(roving{i},[floor(min_rove)+0.5:ceil(max_rove)-0.5]); % center bins
            if do_plot
                figure
                hist(roving{i},[floor(min_rove)+0.5:ceil(max_rove)-0.5]); % center bins
                ylabel('Frequency')
                xlabel('Roving range [dB]')
            end
%             x_dB = x_dB + 0.5;
            
            switch i
                case 1
                    stRove.roveACE = hist(roving{i},vls_rove);
                    stRove.roveACEScore(dist2minrove:dist2minrove+length(Avg_ans(i,:))-1)=Avg_ans(i,:);
                case 2
                    stRove.roveF0m = hist(roving{i},vls_rove);
                    stRove.roveF0mScore(dist2minrove:dist2minrove+length(Avg_ans(i,:))-1)=Avg_ans(i,:);
            end
        end
    end
end

if(do_plot)
    
    marker_size = 15;
    fontsize    = 16;
    
    figure;
    plot(race, 'MarkerSize', marker_size,'Marker', 's', 'MarkerFaceColor','red', 'LineStyle', '--'); hold on, grid on;
    
    if nargin >= 2
        plot(rf0m, 'MarkerSize', marker_size,'Marker', 's', 'MarkerFaceColor','green', 'LineStyle', '--');
    end

    % chance and significance level
    line([1 4], [50 50], 'LineStyle', '--', 'Color', 'green');
    slev = get_significance_level(TimesPresented, 2);
    line([1 4], [slev slev], 'LineStyle', '--', 'Color', 'red');

    set(gca, 'XTick', 1:4);
    set(gca,'XTickLabel', xticks);
    set(gca, 'FontSize', fontsize);
    ylim([0 110]);
    xlim([0.9 length( xticks )+0.1]);
    xlabel('Reference Tones','FontSize',fontsize );
    ylabel('Percent correct');
    title(['PR results; ref ' ref_note ' / each trial presented ' num2str(TimesPresented) ' times'], 'FontSize', fontsize);
    legend( txtLegend, 'Location', 'NorthEastOutside');
    
    disp(['Results strat1:[%]' Num2Str(race)])
    if nargout == 2
        disp(['Results strat2:[%]' Num2Str(rf0m)])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Roving analysis:
    try
        if length(roving) ~= 0

            idx = find(isnan(Avg_ans));
            Avg_ans(idx) = 0;

            figure
            numRove = strtrim(cellstr(num2str([Avg_num(1,:)'],'(%d)')));

            disp(['Times per rove:   ' Num2Str(Avg_num)])
            disp(['Roving influence: ' Num2Str(Avg_ans)])
            disp(['Mean influence:   ' Num2Str(sum(Avg_ans(1,:).*Avg_num(1,:))/sum(Avg_num(1,:)))])

            Avg_ans(idx) = NaN;       

            plot(x_dB, Avg_ans(1,:),'bx','LineWidth',5), hold on
            text(x_dB,Avg_ans(1,:),numRove,'VerticalAlignment','bottom','FontSize',fontsize-2);
            try
                plot(x_dB, Avg_ans(2,:),'ro','LineWidth',5)
                numRove = strtrim(cellstr(num2str([Avg_num(2,:)'],'(%d)')));
                text(x_dB,Avg_ans(2,:),numRove,'VerticalAlignment','bottom','FontSize',fontsize-2);
                legend('ACE', 'F0mod')
            catch
                legend('Current strategy')
            end
            title(['PR (Roving effect); ref ' ref_note ' / each trial presented ' num2str(TimesPresented) ' times'], 'FontSize', fontsize);
            ylabel('Average Correct/Incorrect')
            xlabel('\Delta dB between test and reference (+\Delta dB means test tone louder)')
            ylim([-0.1 1.1])
            grid on
        end
    catch
        disp([mfilename '.m: The APEX result file being read seems to correspond to an old version of PR, no roving is available for plotting'])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end