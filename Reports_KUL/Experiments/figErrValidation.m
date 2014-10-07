function [error,ValuesSim] = figErrValidation(filename, prefix, suffix, dirExp, dirExpSuf, info)
% function [error,ValuesSim] = figErrValidation(filename, prefix, suffix, dirExp, dirExpSuf, info)
%
% prefix - name of the clean audio file (audio files in quiet, without extension), 
% suffix - '-Errors.mat'
% dirExp - main directory where the variables are located
% dirExpSuf is the subdirectory where the variables are located
%
% Example:
% 
% adapted from mergeStruct.m
% Programmed by Alejandro Osses, ExpORL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ValuesSim = [];
error = [];

% List_of_files = {   'confidenceInterval.m', ...
%                     'setupPlotConf.m'}; 
% 
% [AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files,'alejandro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bSaveFigures = 1;

% if nargin < 5
%     info.plotMain = 1;
%     info.bPlotSub = 0;
% end
% 
% if info.bPlotSub == 0
%     info.plotMain = 1;
% end
% 
% if ~isfield(info, 'bPlotSub')
%     info.bPlotSub = 0;
%     info.plotMain = 1;
% end
% 
% if info.bPlotSub == 1
%     info.plotMain = 0;
% end
% 
% if nargin == 0
%     prefix = {  'Choice',...
%                 'wdz6', ...
%                 'Sweep-log-44100-up-to-1kHz-20dBFS', ...
%                 'msine-Gsh2-44100-20dBFS', ...
%                 'msine-C4-44100-20dBFS'};
% else
%     if length(prefix) == 0
%         prefix = {  'Choice',...
%                 'wdz6', ...
%                 'Sweep-log-44100-up-to-1kHz-20dBFS', ...
%                 'msine-Gsh2-44100-20dBFS', ...
%                 'msine-C4-44100-20dBFS'};
%     end
% end
% 
% if length(prefix) == 1
%     stPlot.TitleSuffix = prefix{1};
% end
% 
% if nargin < 2
%     suffix = '-Errors.mat';
% end
% 
% if nargin < 3
%     dirExp = '/home/alejandro/Documenten/LaTeX_Docs/predoc_plan_wekelijkse_update/wu026_2013_06_05/Variables/';
% end
% 
% try
%     dirExp = info.ErrorMoveTo;
% catch
% 
%     try
%         dirExp = [dirExp dirExpSuf{1}];
%     end
%     try
%         dirExp = [dirExp dirExpSuf{2}];
%     end
%     try
%         dirExp = [dirExp dirExpSuf{3}];
%     end
%         
% end

ref = {'m05' 'p00' 'p05' 'p10' 'p20' 'Q'; ...
        '-5'  '0'  '+5'  '+10' '+20' 'Q';
        1     2     3     4     5     6};
try % loading errValidation variable
    load(filename)
catch
    load([info.root_location filename])
end

for k = 1:length(errValidation) % One per each audio file
%     
%     NameCompare = deleteExtension( prefix{k} );
%     
%     
%     for i = 1:length(list) % Load every 'error-file'
%         name = list(i).name;
        
        
%         load([dirExp name]);
        start = 1; %length(name)-2-length(suffix);
%         ID = name(start:start+2);
% 
%         idx = 0;
        for m = 1:5
%             if strcmp(ID,ref{1,m})
%                 idx = ref{3,m};
%             end
        end
%         if idx == 0
%             idx = 6;
%         end
% 
        errValidation(1,1).Sim(1,1).vErr
        errValidation(1,1).Sim(1,1).uvErr
        % other audio
        errValidation(2,1).Sim(1,1).vErr
        errValidation(2,1).Sim(1,1).uvErr
        
        error.MatrixSim.vErr(k,idx)     = errTot.Sim.vErr;
        error.MatrixSim.uvErr(k,idx)    = errTot.Sim.uvErr;
        error.MatrixSim.uvErrTotal(k,idx) = errTot.Sim.uvErrTotal;
        error.MatrixSim.gErr(k,idx)     = errTot.Sim.gErr;
        error.MatrixSim.f0Dev(k,idx)    = errTot.Sim.f0Dev;
%     end
%     
end
% 
% %%
% nCols       = 2;
% xTick       = 1:size(ref,2); % Defines the x-axis on plot, size: 1x4 
% xTickLabel  = {ref{2,:}};
% xLim        = [0.8      6.2]; % Manual setup
% offset_y    = 0;
% yLim(1,:)   = [-50-offset_y      100+offset_y]; % between -60 and 60
% yTick{1}    = [yLim(1,1)+offset_y+10:10:yLim(1,2)-offset_y-10];
% yLim(2,:)   = [-28-offset_y      2*28+offset_y];
% yTick{2}    = [yLim(2,1)+offset_y+10:10:yLim(2,2)+offset_y-10];
% yTickLabel  = {abs(yTick{1})};
% XLabel      = 'SNR (dB)';
% YLabel      = {'uvErr (%)            vErr (%)',...
%                'f0Dev (%)           gErr (%)'};
% SeriesLabel = {'NMT', 'Sim'};
% 
% stPlot.SeriesLabel = SeriesLabel;
% stPlot.xTick    = xTick;
% stPlot.yTick    = yTick;
% stPlot.xTickLabel = xTickLabel;
% stPlot.yTickLabel = yTickLabel;
% stPlot.xLim     = xLim;
% stPlot.yLim     = yLim;
% stPlot.TitleHead = '';
% stPlot.XLabel = XLabel;
% stPlot.YLabel = YLabel;
% stPlot.ReverseData = 1;
% stPlot.nCols = nCols;
% stPlot.sameYAxis = 0;
% stPlot.LocationLegend = 'NorthEast';
% stPlot.figPos  = [0 0 750 768/2];
% % stPlot.figPos  = [0 0 1000 768/2];
% 
% i_min = 1; %255;
% i_max = length(prefix); %264;
% 
% ac = PlotMeans( [error.MatrixNM.vErr(i_min:i_max,:)+offset_y       error.MatrixNM.gErr(i_min:i_max,:)+offset_y ], ...
%                 [error.MatrixSim.vErr(i_min:i_max,:)+offset_y      error.MatrixSim.gErr(i_min:i_max,:)+offset_y], ...
%                  stPlot, ...
%               [-1*error.MatrixNM.uvErr(i_min:i_max,:)-offset_y   -1*error.MatrixNM.f0Dev(i_min:i_max,:)-offset_y], ...
%               [-1*error.MatrixSim.uvErr(i_min:i_max,:)-offset_y  -1*error.MatrixSim.f0Dev(i_min:i_max,:)-offset_y]); % "Plotted downwards"
%           
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This was used for IWT Meeting:
% %
% % ac = PlotMeans( [error.MatrixNM.vErr(i_min:i_max,:)+offset_y        error.MatrixNM.f0Dev(i_min:i_max,:)+offset_y ], ... 
% %                 [error.MatrixNM.vErr(i_min:i_max,:)-offset_y        error.MatrixSim.f0Dev(i_min:i_max,:)-offset_y], ...
% %                  stPlot, ...
% %               [-1*error.MatrixNM.uvErr(i_min:i_max,:)-offset_y   -50-0*error.MatrixNM.f0Dev(i_min:i_max,:)-offset_y], ...
% %               [-1*error.MatrixSim.uvErr(i_min:i_max,:)-offset_y  -50-0*error.MatrixSim.f0Dev(i_min:i_max,:)-offset_y]);
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if bSaveFigures == 1
%     
%     if info.plotMain % Plot with deviations
%             
%         FiguurNaam = '-validation-F0';
% 
%         try
%             FiguurMap   = [dirExp, num2str(dirExpSuf)];
%             saveas(ac.handle,[FiguurMap, FiguurNaam],'epsc');
%         catch
%             if ~isfield(info,'ErrorMoveTo')
%                 FiguurMap = [dirExp dirExpSuf{1}];
%                 saveas(ac.handle,[FiguurMap, FiguurNaam],'epsc') % Algemeen-validation
%             else
%                 FiguurMap = info.ErrorMoveTo;
%                 FiguurNaam = [dirExpSuf{1}(1:end-1) '-validation-F0'];
%                 saveas(ac.handle,[FiguurMap FiguurNaam],'epsc');
%             end
% 
%         end
%             
%     end
%     
% end    
%     
% XLabel      = 'SNR (dB)';
% YLabel      = {'uvErrTotal (%)       vErr (%)',...
%                'f0Dev (%)           gErr (%)'};
% stPlot.YLabel = YLabel;
% 
% if info.bPlotSub == 1
%     if length(prefix)==1
%         
%         if strcmp(prefix{1}(end-2:end),'wav')
%             FiguurNaam = [prefix{1}(1:end-4) '-validation-F0'];
%         else
%             FiguurNaam = [prefix{1} '-validation-F0'];
%         end
%         FiguurMap = dirExp;
%         if ~isfield(info,'MoveTo')
%             saveas(ac.handle,[FiguurMap,  FiguurNaam],'epsc');
%         else
%             
%             FiguurMap = info.MoveTo;
%             saveas(ac.handle,[FiguurMap FiguurNaam],'epsc');
%         end
%     end
%     
% end   
% 
% display(['Plot saved as: ' FiguurNaam ' in: ' FiguurMap])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if CountAddedPaths ~= 0
%     for i = 1:CountAddedPaths
%         rmpath( AddedPaths{CountAddedPaths} )
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function string_new = deleteExtension(string)
% 
% if strcmp( string(end-2:end),'wav')
%     string_new = string(1:end-4);
% else
%     string_new = string;
% end