% PR_Read_apexresult_Test_Inicial.m
% 
% Reads and processes APEX Result Files (*.apr).
% The experiment files should have the format '<xsltscript>apexresult.xsl</xsltscript>'
% for a PR Test

clear, close all, clc

% Check of dependencies:
% List_of_files = {   'get_mci_contours.m'};
% [AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files);

display(['Starts: ',mfilename])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fontsize    = 10;
DirResults      = '/home/alejandro/Documenten/Meas/Meas/Experiments/Results_XML/AOV_2013_Tests/';

PR_Group1_dir   = { [DirResults 'PR_Freq_131_edit_RZ-results.apr'  ], ...
                    [DirResults 'PR_Freq_131_edit_RZ-results.apr'  ]};
%                   [DirResults 'MCI_AC_104_groep_1-results-1.apr']};
% MCI_Group2_dir = {[DirResults 'MCI_AC_104_groep_2-results.apr'  ], ...
%                   [DirResults 'MCI_AC_104_groep_2-results-1.apr']};
% MCI_Group3_dir = {[DirResults 'MCI_AC_104_groep_3-results.apr'  ], ...
%                   [DirResults 'MCI_AC_104_groep_3-results-1.apr']};
% MCI_Group4_dir = {[DirResults 'MCI_AC_104_groep_4-results.apr'  ], ...
%                   [DirResults 'MCI_AC_104_groep_4-results-1.apr']};
% 
% % Group 1: 1 Semitone %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
p           = PR_Read_apexresult( PR_Group1_dir , 13);
% means_f0m1  = p{1}.Mean;
%  
% % Group 2: 2 Semitones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% p           = MCI_Read_apexresult( MCI_Group2_dir );
% means_f0m2  = p{1}.Mean;
% 
% % Group 3: 3 Semitones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% p           = MCI_Read_apexresult( MCI_Group3_dir );
% means_f0m3  = p{1}.Mean;
% 
% % Group 4: 4 Semitones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% p           = MCI_Read_apexresult( MCI_Group4_dir );
% means_f0m4  = p{1}.Mean;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot:
% 
% figure
% i = 1;
% cmap = [0.3     0.3     0.3     ; ... % Color map for figure
%         0.7460  0.7460  0.7460  ; ...
%         1       1       1       ; ...
%         0.6     0.6     0.6];
%     
% % bar([means_ace means_f0m]);hold on;
% bar([means_f0m1 means_f0m2 means_f0m3 means_f0m4]);hold on;
% colormap(cmap);
% 
% set( gca, 'XTick', 1:nContours );
% set( gca,'XTickLabel', contours );
% set( gca, 'FontSize', fontsize );
%  
% % if (i==1)
% %      title( [num2str(i) ' Semitone'], 'FontSize', fontsize );
% % else
% %  	title( [num2str(i) ' Semitones'], 'FontSize', fontsize );
% % end
% % legend({'ACE', 'F0mod'}, 'Location', 'NorthEastOutside');
% 
% legend({'Group 1', 'Group 2','Group 3', 'Group 4'}, 'Location', 'NorthEastOutside');
%  
% % % slev = get_significance_level( nSubjects*repetitions, 9 );
% % % line( [0 nContours+1], [slev slev], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 2, 'LineStyle', '--' );
% ylim([0 110]);
% ylabel('% Correct')

display(['Ends: ',mfilename])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end