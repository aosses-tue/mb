function experiment_report_20131104_roving
% function experiment_report_20131104_roving
%
% Programmed by Alejandro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

main_folder = '/home/alejandro/Documenten/Meas/Meas/Experiments/Experiment_reports/20131104-Roving/';
dest_folder = [main_folder 'Figures/'];

Times_presented = [ 2 5 9 2 11 2 6 7; ...
                    2 4 8 5  8 9 8 0; ...
                    3 6 2 6  7 8 8 4; ...
                    3 3 5 5  9 8 5 6];
Centre_values   = -3.5:1:3.5;

bSave = 1;

figure
bar(Centre_values, sum(Times_presented));

ylabel('Frequency')
xlabel('Roving range [dB]')

for i = 1:length(Centre_values)
    Label{i} = [num2str(floor(Centre_values(i))) ' - ' num2str(ceil(Centre_values(i)))];
end

set(gca,'XTick',Centre_values)
set(gca,'XTickLabel',Label)

if bSave == 1
    saveas(gca, [dest_folder 'hist_experiments'],'epsc')
end

quick_PR([main_folder 'PR_Ref_131_UW_ACE-LB-lot-repetitions-test.apr'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end