function [error,ValuesSim, ValuesNM, h] = figF0Measures(filenames, directory, dirExpSuf, info)
% function [error,ValuesSim, ValuesNM, h] = figF0Measures(filenames, directory, dirExpSuf, info)
%
% filenames - if it is a cell array: name of the clean audio file (audio files in quiet, without extension), 
%           - if it is a struct array: errors for each audio file 
% suffix - '-Errors.mat'
% directory - main directory where the error variables are located
% dirExpSuf is the subdirectory where the variables are located
%
%   h - corresponds to the handles of the Figures in the following order:
%       vErr, uvErr, gErr, Weighted Measure, flat measure
%
% % Example:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ValuesSim = [];
ValuesNM = [];
h = [];
List_of_files = {   'confidenceInterval.m', ...
                    'setupPlotConf.m'}; 

[AddedPaths, CountAddedPaths] = Check_dependencies_ExpORL(List_of_files,'alejandro');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % if nargin < 5
% %     info.plotMain = 1;
% %     info.bPlotSub = 0;
% % end
% % 
% % if info.bPlotSub == 0
% %     info.plotMain = 1;
% % end
% % 
% % if ~isfield(info, 'bPlotSub')
% %     info.bPlotSub = 0;
% %     info.plotMain = 1;
% % end
% % 
% % if info.bPlotSub == 1
% %     info.plotMain = 0;
% % end
% % 
% % if nargin == 0
% %     filenames = {  'Choice',...
% %                 'wdz6', ...
% %                 'Sweep-log-44100-up-to-1kHz-20dBFS', ...
% %                 'msine-Gsh2-44100-20dBFS', ...
% %                 'msine-C4-44100-20dBFS'};
% % else
% %     if length(filenames) == 0
% %         filenames = {  'Choice',...
% %                 'wdz6', ...
% %                 'Sweep-log-44100-up-to-1kHz-20dBFS', ...
% %                 'msine-Gsh2-44100-20dBFS', ...
% %                 'msine-C4-44100-20dBFS'};
% %     end
% % end
% % 
% % if length(filenames) == 1
% %     stPlot.TitleSuffix = filenames{1};
% % end
% % 

if nargin < 3
    if isunix
        directory = '/home/alejandro/Documenten/LaTeX_Docs/predoc_plan_wekelijkse_update/wu026_2013_06_05/Variables/';
    else
        % directory = '/home/alejandro/Documenten/LaTeX_Docs/predoc_plan_wekelijkse_update/wu026_2013_06_05/Variables/';
    end
end

% % try
% %     directory = info.ErrorMoveTo;
% % %     directory = [directory num2str(dirExpSuf)];
% % catch
% % 
% %     try
% %         directory = [directory dirExpSuf{1}];
% %     end
% %     try
% %         directory = [directory dirExpSuf{2}];
% %     end
% %     try
% %         directory = [directory dirExpSuf{3}];
% %     end
% %         
% % end
% % 
ref = {'m05' 'p00' 'p05' 'p10' 'p20' 'Q'; ...
        '-5'  '0'  '+5'  '+10' '+20' 'Q'; ...
        1     2     3     4     5     6};
% %     
% % bSaveFigures = 1;

% The following is going to be done in case filenames is a cell array:
if iscell(filenames)
    for k = 1:length(filenames)
        
        suffix = '-Errors.mat';
        
        NameCompare = filenames{k};
        
        if strcmp( NameCompare(end-2:end),'wav')
            NameCompare = NameCompare(1:end-4);
        end
        
        list = dir([directory NameCompare '*' suffix]);
    % %     for i = 1:length(list)
    % %         name = list(i).name;
    % %         load([directory name]);
    % %         start = length(name)-2-length(suffix);
    % %         ID = name(start:start+2);
    
    % %         idx = 0;
    % %         for m = 1:5
    % %             if strcmp(ID,ref{1,m})
    % %                 idx = ref{3,m};
    % %             end
    % %         end
    % %         if idx == 0
    % %             idx = 6;
    % %         end
    
    % %         error.MatrixSim.vErr(k,idx)     = errTot.Sim.vErr;
    % %         error.MatrixSim.uvErr(k,idx)    = errTot.Sim.uvErr;
    % %         error.MatrixSim.uvErrTotal(k,idx) = errTot.Sim.uvErrTotal;
    % %         error.MatrixSim.gErr(k,idx)     = errTot.Sim.gErr;
    % %         error.MatrixSim.f0Dev(k,idx)    = errTot.Sim.f0Dev;
            
    % %         error.MatrixNMT.vErr(k,idx)      = errTot.NMT.vErr;
    % %         error.MatrixNMT.uvErr(k,idx)     = errTot.NMT.uvErr;
    % %         error.MatrixNMT.uvErrTotal(k,idx) = errTot.NMT.uvErrTotal;
    % %         error.MatrixNMT.gErr(k,idx)      = errTot.NMT.gErr;
    % %         error.MatrixNMT.f0Dev(k,idx)     = errTot.NMT.f0Dev;
    % %     end
        
    end
    
    min_ylim    = -50;
    max_ylim    = 100;

end

if isstruct(filenames) % in this case we just need to reassign errors according the required nomenclature
    total_time = 0;
    [errResult errSim errNMT err2Plot] = F0mod_Physical_validation_analysis_errors(filenames);
    for k = 1:size(filenames,1)
        for idx = 1:size(filenames(1,1).Sim,2)
            errTot.Sim = filenames(k,1).Sim(1,idx);
            
            % vErr, uvErr, gErr to Seconds
            error.MatrixSim.vErr(k,idx)     = errTot.Sim.vErr /100 * errTot.Sim.Total_voiced;
            error.MatrixSim.uvErr(k,idx)    = errTot.Sim.uvErr/100 * errTot.Sim.Total_unvoiced;
            error.MatrixSim.gErr(k,idx)     = errTot.Sim.gErr /100 * errTot.Sim.Total_voiced;
            error.MatrixSim.f0Dev(k,idx)    = errTot.Sim.f0Dev;
            
            % vErr, uvErr, gErr to Seconds
            errTot.NMT = filenames(k,1).NMT(1,idx);
            error.MatrixNMT.vErr(k,idx)     = errTot.NMT.vErr /100 * errTot.NMT.Total_voiced;
            error.MatrixNMT.uvErr(k,idx)    = errTot.NMT.uvErr/100 * errTot.NMT.Total_unvoiced;
            error.MatrixNMT.gErr(k,idx)     = errTot.NMT.gErr /100 * errTot.Sim.Total_voiced;
            error.MatrixNMT.f0Dev(k,idx)    = errTot.NMT.f0Dev;
            
            if isnan(error.MatrixSim.vErr(k,idx))
                error.MatrixSim.vErr(k,idx) = 0;
            end
            if isnan(error.MatrixSim.uvErr(k,idx))
                error.MatrixSim.uvErr(k,idx) = 0;
            end
            if isnan(error.MatrixSim.gErr(k,idx))
                error.MatrixSim.gErr(k,idx) = 0;
            end
            if isnan(error.MatrixNMT.vErr(k,idx))
                error.MatrixNMT.vErr(k,idx) = 0;
            end
            if isnan(error.MatrixNMT.uvErr(k,idx))
                error.MatrixNMT.uvErr(k,idx) = 0;
            end
            if isnan(error.MatrixNMT.gErr(k,idx))
                error.MatrixNMT.gErr(k,idx) = 0;
            end
           
        end
        total_time = total_time + errTot.Sim.TotalTime;
    end
    avg_total_time = total_time / size(filenames,1);
    % vErr, uvErr, gErr to Percentage, considering average total time.
    error.MatrixSim.vErr(k,:)   = error.MatrixSim.vErr(k,:)/avg_total_time * 100;
    error.MatrixSim.uvErr(k,:)  = error.MatrixSim.uvErr(k,:)/avg_total_time * 100;
    error.MatrixSim.gErr(k,:)   = error.MatrixSim.gErr(k,:)/avg_total_time * 100;
    
    error.MatrixNMT.vErr(k,:)   = error.MatrixNMT.vErr(k,:)/avg_total_time * 100;
    error.MatrixNMT.uvErr(k,:)  = error.MatrixNMT.uvErr(k,:)/avg_total_time * 100;
    error.MatrixNMT.gErr(k,:)   = error.MatrixNMT.gErr(k,:)/avg_total_time * 100;
    
    min_ylim    = -5;
    max_ylim    = 16;
end

%%
nCols       = 2;
xTick       = 1:size(ref,2); % Defines the x-axis on plot, size: 1x4 
xTickLabel  = {ref{2,:}};
xLim        = [0.8      6.2]; % Manual setup
offset_y    = 0;
yLim(1,:)   = [min_ylim-offset_y      max_ylim+offset_y]; % between -60 and 60
yTick{1}    = [yLim(1,1)+offset_y+10:10:yLim(1,2)-offset_y-10];
yLim(2,:)   = [min_ylim-offset_y      max_ylim+offset_y]; %[-28-offset_y      2*28+offset_y];
yTick{2}    = [yLim(2,1)+offset_y+10:10:yLim(2,2)+offset_y-10];
yTickLabel  = {abs(yTick{1})};
SeriesLabel = {'NMT', 'Sim'};

stPlot.SeriesLabel = SeriesLabel;
stPlot.xTick    = xTick;
stPlot.yTick    = yTick;
stPlot.xTickLabel = xTickLabel;
% stPlot.yTickLabel = yTickLabel;
stPlot.xLim     = xLim;
stPlot.yLim     = yLim;
stPlot.TitleHead = '';
stPlot.XLabel = 'SNR (dB)';
stPlot.YLabel = {'uvErr (%)            vErr (%)', 'f0Dev (%)           gErr (%)'};
stPlot.ReverseData = 1;
stPlot.nCols = nCols;
stPlot.sameYAxis = 0;
stPlot.LocationLegend = 'NorthEast';
stPlot.figPos  = [0 0 750 768/2];
stPlot.nPlots = 2;
% stPlot.figPos  = [0 0 1000 768/2];

if iscell(filenames)
    
    i_min = 1; %255;
    i_max = length(filenames); %264;

    ac = PlotMeans( [error.MatrixNMT.vErr(i_min:i_max,:)+offset_y      error.MatrixNMT.gErr(i_min:i_max,:)+offset_y ], ...
                    [error.MatrixSim.vErr(i_min:i_max,:)+offset_y      error.MatrixSim.gErr(i_min:i_max,:)+offset_y ], ...
                     stPlot, ...
                  [-1*error.MatrixNMT.uvErr(i_min:i_max,:)-offset_y  -1*error.MatrixNMT.f0Dev(i_min:i_max,:)-offset_y], ...
                  [-1*error.MatrixSim.uvErr(i_min:i_max,:)-offset_y  -1*error.MatrixSim.f0Dev(i_min:i_max,:)-offset_y]); % "Plotted downwards"
else
    stPlot.nPlots   = 1;
    stPlot.nCols    = 1;
    stPlot.YLabel   = {'vErr (%)'};
    stPlot.yTick    = [yLim(1,1)+offset_y+2:2:yLim(1,2)-offset_y];
    PlotMeans( error.MatrixNMT.vErr, error.MatrixSim.vErr, stPlot );
    h(end+1) = gcf;
    
    title(['F0 performance measures considering ' num2str(size(error.MatrixNMT.vErr,1)) ' audio files'])
        
    stPlot.YLabel   = {'uvErr (%)'};
    PlotMeans( error.MatrixNMT.uvErr, error.MatrixSim.uvErr, stPlot );
    title(['F0 performance measures considering ' num2str(size(error.MatrixNMT.vErr,1)) ' audio files'])
    h(end+1) = gcf;
    
    stPlot.YLabel   = {'gErr (%)'};
    stPlot.yLim     = [-1 1];
    PlotMeans( error.MatrixNMT.gErr, error.MatrixSim.gErr, stPlot );
    title(['F0 performance measures considering ' num2str(size(error.MatrixNMT.vErr,1)) ' audio files'])
    h(end+1) = gcf;
    
    Weight = [.25 .50 .25]';

    stPlot.YLabel   = {'Weighted F0 errors (%)'};
    

    stPlot.yLim     = [-10 34];
    stPlot.yTick    = [stPlot.yLim(1,1)+offset_y+2:( stPlot.yLim(1,2)-stPlot.yLim(1,1) )/10:stPlot.yLim(1,2)-offset_y];
    PlotMeans( err2Plot.WeightedNMT, err2Plot.WeightedSim, stPlot )
    title(['F0 performance measures considering ' num2str(size(error.MatrixNMT.vErr,1)) ' audio files'])
    h(end+1) = gcf;
    
    stPlot.YLabel   = {'F0 errors (%) (no weighting)'};
    PlotMeans( err2Plot.FlatNMT, err2Plot.FlatSim, stPlot );
    title(['F0 performance measures considering ' num2str(size(error.MatrixNMT.vErr,1)) ' audio files'])
    h(end+1) = gcf;
end
% %           
% % if bSaveFigures == 1
% %     
% %     if info.plotMain % Plot with deviations
% %             
% %         FiguurNaam = '-validation-F0';
% % 
% %         try
% %             FiguurMap   = [directory, num2str(dirExpSuf)];
% %             saveas(ac.handle,[FiguurMap, FiguurNaam],'epsc');
% %         catch
% %             if ~isfield(info,'ErrorMoveTo')
% %                 FiguurMap = [directory dirExpSuf{1}];
% %                 saveas(ac.handle,[FiguurMap, FiguurNaam],'epsc') % Algemeen-validation
% %             else
% %                 FiguurMap = info.ErrorMoveTo;
% %                 FiguurNaam = [dirExpSuf{1}(1:end-1) '-validation-F0'];
% %                 saveas(ac.handle,[FiguurMap FiguurNaam],'epsc');
% %             end
% % 
% %         end
% %             
% %     end
% %     
% % end    
% %     
% % XLabel      = 'SNR (dB)';
% % YLabel      = {'uvErrTotal (%)       vErr (%)',...
% %                'f0Dev (%)           gErr (%)'};
% % stPlot.YLabel = YLabel;
% % 
% % if info.bPlotSub == 1
% %     if length(filenames)==1
% %         
% %         if strcmp(filenames{1}(end-2:end),'wav')
% %             FiguurNaam = [filenames{1}(1:end-4) '-validation-F0'];
% %         else
% %             FiguurNaam = [filenames{1} '-validation-F0'];
% %         end
% %         FiguurMap = directory;
% %         if ~isfield(info,'MoveTo')
% %             saveas(ac.handle,[FiguurMap,  FiguurNaam],'epsc');
% %         else
% %             
% %             FiguurMap = info.MoveTo;
% %             saveas(ac.handle,[FiguurMap FiguurNaam],'epsc');
% %         end
% %     end
% %     
% % end   
% % 
% % display(['Plot saved as: ' FiguurNaam ' in: ' FiguurMap])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CountAddedPaths ~= 0
    for i = 1:CountAddedPaths
        rmpath( AddedPaths{CountAddedPaths} )
    end
end

end