function [race, rf0m, labelContours] = plot_MC_results( filename, InitialSubject, ref_note, subject, do_plot )
                                     % ( root_note, interval, splits, subject, repetitions, plot_details, do_plot)

% function [race, rf0m, labelContours] = plot_MC_results( root_note, interval, splits, subject, repetitions, plot_details, do_plot)
% Outputs:
%   res_ace_all
%   res_f0m_all
%   conf_ace
%   conf_f0m
%
% Programmed by Alejandro, adapted from a Matthias' script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
	return;
end

if nargin ~= 1
    result_files = {[filename 'ACE-' InitialSubject '.apr'], ...
                    [filename 'F0m-' InitialSubject '.apr']};
	txtLegend = {'ACE', 'F0mod'};
else
    result_files = { filename };
    do_plot = 1;
    txtLegend   = {'Instant result'};
    subject     = ': Current subject';
end
            
nContours       = 9;
splits          = 1; % Non-splitted % To delete for future releases
interval        = 1;

try
    [h, xx, xx, xx, TimesPresented, nTotalTrials] = eval_results_MC(result_files, 0, txtLegend);
catch
    result_files = {[filename(1:end-1) '_ACE-ACE-' InitialSubject '.apr'], ... % Loudness balanced experiments
                    [filename(1:end-1) '_F0m-F0m-' InitialSubject '.apr']};
	[h, xx, xx, xx, TimesPresented, nTotalTrials] = eval_results_MC(result_files, 0, txtLegend);
end

nRoving         = TimesPresented;
nRepetitions    = 1;

hace = h{1}{1};

if nargout >= 2
    hf0m = h{2}{1};
end

race = zeros(1,9);

if nargout
    rf0m = zeros(1,9);
end

xticks = get_mci_contours(1); % Labels for each Contour

for i=1:nContours
    
    for j=1:size(hace, 1)
        if( strcmp(hace{j,1}(1:end-1), xticks{i} )) % hace{j,1}(1:end-1)
            race(i) = hace{j, 2}/TimesPresented*100; % percentage
        end
    end
    
    if nargout >= 2
        for j=1:size(hf0m, 1)
            if( strcmp(hf0m{j,1}(1:end-1), xticks{i} ))
                rf0m(i) = hf0m{j, 2}/TimesPresented*100; % percentage
            end
        end
    end

end

if(do_plot)
    
    marker_size = 15;
    fontsize    = 16;

    figure;
    plot(   1:length(race), race, 'MarkerSize', marker_size,'Marker', 's', 'MarkerFaceColor','red', 'LineStyle', '--' ); hold on, grid on;
    
    if nargout >=2
        plot(   1:length(rf0m), rf0m );
    end
    
    set(gca, 'XTick', 1:nContours);
    set(gca,'XTickLabel',xticks)
    
    set(gca, 'FontSize', fontsize);
    ylim([0 110]);
    title([subject ' MCI Int ' interval ' sessions ' num2str(splits)]);
    
end

labelContours = xticks;

end
