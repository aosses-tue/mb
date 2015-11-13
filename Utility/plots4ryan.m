function plots4ryan
% function plots4ryan
%
% 1. Description:
%
% 2. Stand-alone example:
%
% 3. Additional info:
%       Tested cross-platform: No
%
% Programmed by Alejandro Osses, HTI, TU/e, the Netherlands, 2014-2015
% Created on    : 11/11/2015
% Last update on: 11/11/2015 
% Last use on   : 11/11/2015 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
pathfigs = 'D:\Downloads\';

figs = {'Precedence_02.fig', ...
        'Precedence_03.fig', ...
        'MOC_02.fig', ...
        'MOC_03.fig'};
h = [];
for i = 1:length(figs)
    open([pathfigs figs{i}]);
    h(end+1) = gcf;
    
    if i == 4
        ylim([-4 5]);
    end
end

h2save(1) = Figure2paperfigure(h(1:2),2,1);
h2save(2) = Figure2paperfigure(h(3:4),2,1);

outputfilename = [pathfigs 'merged-fig-precedence'];

Saveas(h2save(1),outputfilename,'epsc'); % eps color
Saveas(h2save(1),outputfilename,'fig'); % eps color

outputfilename = [pathfigs 'merged-fig-MOC'];

Saveas(h2save(2),outputfilename,'epsc'); % eps color
Saveas(h2save(2),outputfilename,'fig'); % eps color

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
