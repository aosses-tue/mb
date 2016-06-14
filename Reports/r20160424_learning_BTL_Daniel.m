function r20160424_learning_BTL_Daniel
% function r20160424_learning_BTL_Daniel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

labels = {'LBJ','HW','CdG','JU','CY','AJF','BB','ET','SL'};
    
scores = [  0 159 163 175 183 179 173 160 142; ...
            75 0 138 164 172 160 156 122 122; ...
            71 96 0 145 157 138 140 122 120; ...
            59 70 89 0 176 115 124 86 61; ...
            51 62 77 58 0 77 95 72 61; ...
            55 74 96 119 157 0 134 92 71; ...
            61 78 94 110 139 100 0 67 48; ...
            74 112 112 148 162 142 167 0 87; ...
            92 112 114 173 173 163 186 147 0];
var2latex(scores);

for i = 1:size(scores,1)
    scores(i,i)=nan;
end

[M_col_over_row S_col_over_row] = Get_mean(scores);
[M_row_over_col S_row_over_col] = Get_mean(scores');
% S_col_over_row = std(scores);
% S_row_over_col = std(scores');


var2latex(round([M_col_over_row; S_col_over_row]*10)/10)

var2latex(round(10*[ M_row_over_col' S_row_over_col'])/10)

Npairs    = 36;
Nsubjects = 234;

Nscores = sum(sum(scores)); % 36 x 234 = 8424

scoresN = scores/Nsubjects;

Pij = 0.02:.05:1;
Dij = il_BTL(Pij);

dirout = 'D:\Documenten-TUe\01-Text\05-Doc-TUe\lx2016-04-21-update\Figures\';
FontSize = 14;
figure;
plot(Pij,Dij,'o--'); grid on
Xlabel('Input score P_i_j',FontSize);
Ylabel('Distance score D_i_j',FontSize);
Title('BTL scale',FontSize)
set(gca,'FontSize',FontSize)
Saveas(gcf,[dirout 'BTL-scale'],'epsc');

figure;
plot(il_BTL(M_row_over_col/Nsubjects),'s');
ha = gca;
set(ha,'XTick',1:length(M_row_over_col))
set(ha,'XTickLabel',labels)
xlim([0.5 length(M_row_over_col)])
grid on;
Ylabel('Distance D_i_j')
Xlabel('Preferred stimulus (labels)')
Title('BTL scale applied to the example',FontSize)
set(gca,'FontSize',FontSize)
Saveas(gcf,[dirout 'BTL-scale-example'],'epsc');


disp('')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = [47 	51 	49 	52 	51 	46 	63 	64; ...
	46 	53 	50 	56 	46 	48 	70 	62; ...
	50 	57 	42 	46 	46 	47 	63 	66; ...
	52 	54 	48 	52 	45 	55 	58 	64; ...
	46 	55 	60 	53 	52 	49 	59 	62; ...
	36 	53 	47 	49 	54 	61 	61 	62; ...
	47 	54 	51 	52 	48 	53 	67 	58; ...
	46 	57 	57 	50 	47 	48 	64 	62; ...
	36 	61 	49 	50 	47 	50 	59 	67; ...
	44 	57 	49 	49 	54 	44 	61 	59];

mean(data)
% Mean 	45.0 	55.2 	50.2 	50.9 	49.0 	50.1 	62.5 	62.6
% Correlation 	r=-0.33 	r=0.49 	r=0.06 	r=-0.39
% Significance 	P=0.35 	P=0.15 	P=0.86 	P=0.27 

xlist = data(:,1);
ylist = data(:,2);
p=polyfit(xlist, ylist, 1);

disp('')
function Dij = il_BTL(Pij)

Dij = log(Pij) - log(1-Pij);